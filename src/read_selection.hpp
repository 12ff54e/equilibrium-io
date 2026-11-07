#pragma once

#include "equilibrium/binary.hpp"
#include "equilibrium/validation.hpp"

#include <stdexcept>

namespace equilibrium::io_detail {

class ReadSelection {
   public:
    explicit ReadSelection(const ReadRequest& request, bool strict = false)
        : request_(request), strict_(strict) {}

    bool wants(DataGroup group) const {
        return request_[group] != Requirement::IGNORED;
    }
    bool accept(DataGroup group, bool present, const std::string& error = {}) {
        auto& status = diagnostics[group];
        if (!wants(group)) return false;
        status.availability = !error.empty() ? Availability::INVALID
                              : present      ? Availability::AVAILABLE
                                             : Availability::ABSENT;
        status.diagnostic = error;
        if (strict_ && status.availability == Availability::INVALID)
            throw std::invalid_argument(error);
        if ((request_[group] == Requirement::REQUIRED &&
             status.availability != Availability::AVAILABLE))
            throw std::runtime_error(
                std::string(group_name(group)) + ": " +
                (error.empty() ? "group is absent" : error));
        return status.availability == Availability::AVAILABLE;
    }
    void unavailable(DataGroup group, const std::string& error) {
        accept(group, false, error);
    }
    void reconcile(Equilibrium& owner) {
        const auto drop = [&](DataGroup group, const std::string& reason) {
            unavailable(group, reason);
            switch (group) {
                case DataGroup::GEOMETRY:
                    owner.families.clear();
                    break;
                case DataGroup::FIELDS:
                    owner.fields.reset();
                    break;
                case DataGroup::GRID:
                    owner.grid.reset();
                    break;
                case DataGroup::INPUT:
                    owner.input.reset();
                    break;
                default:
                    break;
            }
        };
        // Resolve required-group conflicts first, so an optional group's
        // disagreement cannot discard another optional group that agrees with
        // the consumer's required data. Two conflicting optional groups are
        // both omitted; two required groups make the request fail.
        for (bool required_pass : {true, false}) {
            const auto conflict = [&](DataGroup a, DataGroup b, bool differs) {
                if (!differs) return;
                const bool ar = request_[a] == Requirement::REQUIRED;
                const bool br = request_[b] == Requirement::REQUIRED;
                if ((ar || br) != required_pass) return;
                const auto reason = std::string(group_name(a)) + " and " +
                                    group_name(b) + " metadata disagree";
                if (!ar || br) drop(a, reason);
                if (!br) drop(b, reason);
            };
            conflict(DataGroup::GRID, DataGroup::FIELDS,
                     owner.grid && owner.fields &&
                         (owner.grid->ntheta != owner.fields->ntheta ||
                          owner.grid->nzeta != owner.fields->nzeta));
            conflict(DataGroup::GRID, DataGroup::INPUT,
                     owner.grid && owner.input &&
                         *owner.grid != grid_metadata(*owner.input));
            conflict(DataGroup::FIELDS, DataGroup::INPUT,
                     owner.fields && owner.input &&
                         (owner.fields->ntheta != owner.input->ntheta ||
                          owner.fields->nzeta != owner.input->nzeta));
            conflict(DataGroup::GEOMETRY, DataGroup::INPUT,
                     owner.has_geometry() && owner.input &&
                         (owner.lasym() != owner.input->lasym ||
                          static_cast<std::size_t>(owner.input->mpol) *
                                  (static_cast<std::size_t>(owner.input->ntor) +
                                   1) !=
                              static_cast<std::size_t>(owner.mnmax)));
        }
    }
    void finish() {
        for (std::size_t i = 0; i < request_.groups.size(); ++i) {
            const auto group = static_cast<DataGroup>(i);
            if (wants(group) &&
                diagnostics[group].availability == Availability::NOT_REQUESTED)
                accept(group, false);
        }
    }

    ReadDiagnostics diagnostics;

   private:
    const ReadRequest& request_;
    bool strict_;
};

Result<Equilibrium> read_native_binary(const std::string& path,
                                       const ReadRequest& request,
                                       bool strict = false);

inline ReadRequest legacy_request(bool metadata) {
    ReadRequest request;
    request[DataGroup::GEOMETRY] = Requirement::REQUIRED;
    request[DataGroup::FIELDS] = Requirement::OPTIONAL;
    if (metadata) {
        request[DataGroup::GRID] = Requirement::OPTIONAL;
        request[DataGroup::INPUT] = Requirement::OPTIONAL;
        request[DataGroup::PROVENANCE] = Requirement::OPTIONAL;
    }
    return request;
}

inline EquilibriumSnapshot take_snapshot(Equilibrium owner) {
    EquilibriumSnapshot state;
    state.ns = owner.ns;
    state.mnmax = owner.mnmax;
    state.source_version = owner.source.version;
    state.source_checkpoint = owner.source.checkpoint;
    state.families = std::move(owner.families);
    if (owner.fields) {
        state.ntheta = owner.fields->ntheta;
        state.nzeta = owner.fields->nzeta;
        state.half_fields = std::move(owner.fields->half);
        state.full_fields = std::move(owner.fields->full);
    }
    return state;
}

struct ContainerReadSelection {
    ReadSelection selection;
    std::optional<GridMetadata> grid;

    explicit ContainerReadSelection(const ReadRequest& request)
        : selection(request) {}

    Result<Equilibrium> finish(EquilibriumSnapshot state,
                               RunReport report,
                               OutputFormat format) {
        try {
            auto owner = take_equilibrium(std::move(state));
            owner.source.format = format;
            const auto available = [&](DataGroup group) {
                return selection.diagnostics[group].availability ==
                       Availability::AVAILABLE;
            };
            if (!available(DataGroup::GEOMETRY)) owner.families.clear();
            if (!available(DataGroup::FIELDS)) owner.fields.reset();
            if (available(DataGroup::GRID)) owner.grid = grid;
            if (available(DataGroup::INPUT))
                owner.input = std::move(report.input_params);
            if (available(DataGroup::PROVENANCE))
                owner.provenance =
                    std::move(static_cast<RunProvenance&>(report));
            for (auto group :
                 {DataGroup::GEOMETRY, DataGroup::FIELDS, DataGroup::GRID,
                  DataGroup::INPUT, DataGroup::PROVENANCE}) {
                if (!available(group)) continue;
                auto status = validate_group(as_view(owner), group);
                if (selection.accept(group, true, status ? "" : status.error()))
                    continue;
                switch (group) {
                    case DataGroup::GEOMETRY:
                        owner.families.clear();
                        break;
                    case DataGroup::FIELDS:
                        owner.fields.reset();
                        break;
                    case DataGroup::GRID:
                        owner.grid.reset();
                        break;
                    case DataGroup::INPUT:
                        owner.input.reset();
                        break;
                    case DataGroup::PROVENANCE:
                        owner.provenance.reset();
                        break;
                    default:
                        break;
                }
            }
            selection.reconcile(owner);
            selection.finish();
            owner.read_diagnostics = std::move(selection.diagnostics);
            return owner;
        } catch (const std::exception& error) {
            return Result<Equilibrium>(error.what());
        }
    }
};

}  // namespace equilibrium::io_detail
