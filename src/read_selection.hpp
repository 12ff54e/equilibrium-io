#pragma once

#include "equilibrium/binary.hpp"
#include "equilibrium/validation.hpp"

#include <stdexcept>

namespace equilibrium::io_detail {

using FieldShapes = std::array<std::optional<FieldGrid>,
                               static_cast<std::size_t>(Quantity::COUNT)>;

inline void drop_quantity(Equilibrium& owner, Quantity q) {
    if (q == Quantity::GEOMETRY || q == Quantity::SPECTRAL_SHIFT) {
        for (std::size_t i = 0; i < owner.families.size(); ++i)
            if ((i % 3 == 2) == (q == Quantity::SPECTRAL_SHIFT))
                owner.families[i].clear();
    } else if (q == Quantity::GRID)
        owner.grid.reset();
    else if (q == Quantity::INPUT)
        owner.input.reset();
    else if (q == Quantity::PROVENANCE)
        owner.provenance.reset();
    else if (owner.fields) {
        if (q == Quantity::CURRENT_DENSITY)
            for (auto& v : owner.fields->full) v.clear();
        for (std::size_t i = 0; i < owner.fields->half.size(); ++i)
            if (half_field_quantity(i) == q) owner.fields->half[i].clear();
    }
}

class ReadSelection {
   public:
    explicit ReadSelection(const ReadRequest& request) : request_(request) {}

    bool wants(Quantity q) const { return request_[q] != Requirement::IGNORED; }
    bool wants(DataGroup group) const {
        return request_.strongest_requirement(group) != Requirement::IGNORED;
    }
    Requirement requirement(Quantity q) const { return request_[q]; }
    bool available(Quantity q) const {
        return diagnostics[q].availability == Availability::AVAILABLE;
    }
    bool accept(Quantity q, bool present, const std::string& error = {}) {
        if (!wants(q)) return false;
        auto& status = diagnostics[q];
        status.availability = !error.empty() ? Availability::INVALID
                              : present      ? Availability::AVAILABLE
                                             : Availability::ABSENT;
        status.diagnostic = error;
        if (request_[q] == Requirement::REQUIRED && !available(q))
            throw std::runtime_error(
                std::string(quantity_name(q)) + ": " +
                (error.empty() ? "quantity is absent" : error));
        return available(q);
    }
    bool accept(DataGroup group, bool present, const std::string& error = {}) {
        bool any = false;
        for (auto q : group_quantities(group))
            any = accept(q, present, error) || any;
        return any;
    }
    template <class Key>
    void unavailable(Key key, const std::string& error) {
        accept(key, false, error);
    }

    template <class Drop>
    void conflict(std::span<const Quantity> quantities,
                  const std::string& reason,
                  Drop drop) {
        bool required = false;
        for (auto q : quantities)
            required = required || requirement(q) == Requirement::REQUIRED;
        for (auto q : quantities) {
            if (required && requirement(q) == Requirement::REQUIRED) continue;
            unavailable(q, reason);
            drop(q);
        }
        if (required) {
            // If no optional participant could resolve the conflict, all
            // participants are required and the inconsistent read must fail.
            bool all_required = true;
            for (auto q : quantities)
                all_required =
                    all_required && requirement(q) == Requirement::REQUIRED;
            if (all_required) throw std::runtime_error(reason);
        }
    }

    void reconcile(Equilibrium& owner,
                   std::optional<std::reference_wrapper<const FieldShapes>>
                       shapes = std::nullopt) {
        const auto shape = [&](Quantity q) {
            if (shapes)
                return shapes->get()[static_cast<std::size_t>(q)].value_or(
                    FieldGrid{});
            return owner.fields ? FieldGrid{owner.ns, owner.fields->ntheta,
                                            owner.fields->nzeta}
                                : FieldGrid{};
        };
        // Required quantities have priority over optional conflicting data.
        for (bool required_pass : {true, false}) {
            const auto check = [&](Quantity a, Quantity b, bool differs) {
                if (!differs || !available(a) || !available(b)) return;
                const bool required = requirement(a) == Requirement::REQUIRED ||
                                      requirement(b) == Requirement::REQUIRED;
                if (required != required_pass) return;
                const auto reason = std::string(quantity_name(a)) + " and " +
                                    quantity_name(b) + " metadata disagree";
                conflict(std::array{a, b}, reason,
                         [&](Quantity q) { drop_quantity(owner, q); });
            };
            check(Quantity::GRID, Quantity::INPUT,
                  owner.grid && owner.input &&
                      *owner.grid != grid_metadata(*owner.input));
            for (auto q : group_quantities(DataGroup::FIELDS)) {
                const auto a = shape(q);
                for (auto peer : group_quantities(DataGroup::FIELDS)) {
                    const auto b = shape(peer);
                    if (q < peer)
                        check(q, peer,
                              a.ntheta != b.ntheta || a.nzeta != b.nzeta);
                }
                check(Quantity::GRID, q,
                      owner.grid && owner.fields &&
                          (owner.grid->ntheta != a.ntheta ||
                           owner.grid->nzeta != a.nzeta));
                check(q, Quantity::INPUT,
                      owner.fields && owner.input &&
                          (a.ntheta != owner.input->ntheta ||
                           a.nzeta != owner.input->nzeta));
            }
            for (auto q : group_quantities(DataGroup::GEOMETRY))
                check(
                    q, Quantity::INPUT,
                    owner.input &&
                        (owner.lasym() != owner.input->lasym ||
                         static_cast<std::size_t>(owner.input->mpol) *
                                 (static_cast<std::size_t>(owner.input->ntor) +
                                  1) !=
                             static_cast<std::size_t>(owner.mnmax)));
        }
    }

    void finish() {
        for (std::size_t i = 0; i < request_.quantities.size(); ++i) {
            const auto q = static_cast<Quantity>(i);
            if (wants(q) &&
                diagnostics[q].availability == Availability::NOT_REQUESTED)
                accept(q, false);
        }
        for (std::size_t i = 0; i < diagnostics.groups.size(); ++i) {
            auto group = static_cast<DataGroup>(i);
            auto& status = diagnostics[group];
            std::size_t count = 0;
            bool invalid = false, absent = false;
            for (auto q : group_quantities(group)) {
                const auto& member = diagnostics[q];
                count += member.availability == Availability::AVAILABLE;
                invalid =
                    invalid || member.availability == Availability::INVALID;
                absent = absent || member.availability == Availability::ABSENT;
                if (!member.diagnostic.empty()) {
                    if (!status.diagnostic.empty()) status.diagnostic += "; ";
                    status.diagnostic += std::string(quantity_name(q)) + ": " +
                                         member.diagnostic;
                }
            }
            status.availability = count == group_quantities(group).size()
                                      ? Availability::AVAILABLE
                                  : count   ? Availability::PARTIAL
                                  : invalid ? Availability::INVALID
                                  : absent  ? Availability::ABSENT
                                            : Availability::NOT_REQUESTED;
        }
    }

    ReadDiagnostics diagnostics;

   private:
    const ReadRequest& request_;
};

Result<Equilibrium> read_native_binary(const std::string& path,
                                       const ReadRequest& request);

struct ContainerReadSelection {
    ReadSelection selection;
    std::optional<GridMetadata> grid;
    FieldShapes field_shapes;

    explicit ContainerReadSelection(const ReadRequest& request)
        : selection(request) {}

    void record_field_shape(Quantity q, const EquilibriumSnapshot& state) {
        field_shapes[static_cast<std::size_t>(q)] =
            FieldGrid{state.ns, state.ntheta, state.nzeta};
    }

    Result<Equilibrium> finish(EquilibriumSnapshot state,
                               RunReport report,
                               OutputFormat format) {
        try {
            auto owner = take_equilibrium(std::move(state));
            owner.source.format = format;
            if (selection.available(Quantity::GRID)) owner.grid = grid;
            if (selection.available(Quantity::INPUT))
                owner.input = std::move(report.input_params);
            if (selection.available(Quantity::PROVENANCE))
                owner.provenance =
                    std::move(static_cast<RunProvenance&>(report));
            for (std::size_t i = 0;
                 i < static_cast<std::size_t>(Quantity::COUNT); ++i) {
                const auto q = static_cast<Quantity>(i);
                if (!selection.available(q)) {
                    drop_quantity(owner, q);
                    continue;
                }
                if (const auto& shape = field_shapes[i];
                    shape && owner.fields) {
                    owner.fields->ntheta = shape->ntheta;
                    owner.fields->nzeta = shape->nzeta;
                }
                const auto status = validate_quantity(as_view(owner), q);
                if (!selection.accept(q, true, status ? "" : status.error()))
                    drop_quantity(owner, q);
            }
            selection.reconcile(owner, std::cref(field_shapes));
            if (owner.fields)
                for (auto q : group_quantities(DataGroup::FIELDS))
                    if (selection.available(q)) {
                        const auto& shape =
                            field_shapes[static_cast<std::size_t>(q)];
                        owner.fields->ntheta = shape->ntheta;
                        owner.fields->nzeta = shape->nzeta;
                        break;
                    }
            selection.finish();
            owner.read_diagnostics = std::move(selection.diagnostics);
            return owner;
        } catch (const std::exception& error) {
            return Result<Equilibrium>(error.what());
        }
    }
};

}  // namespace equilibrium::io_detail
