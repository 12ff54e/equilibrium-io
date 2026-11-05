import assert from 'node:assert/strict';
import {mkdir, readFile, writeFile} from 'node:fs/promises';
import {join} from 'node:path';
import {createFormatsCpu} from '../js/formats_cpu.mjs';
import {createCpuWasm} from '../js/cpu_wasm.mjs';
import {cpuFiles} from '../js/cpu_files.mjs';

const [fixturePath, productionPath, output] = process.argv.slice(2);
const fixture = await createCpuWasm(await readFile(fixturePath));
const files = cpuFiles(fixture);
const binary = await WebAssembly.compile(await readFile(productionPath));
for (const {module, name} of WebAssembly.Module.imports(binary))
  assert(module === 'wasi_snapshot_preview1' || module === 'env' && name === 'emscripten_notify_memory_growth',
    'Unexpected formats Wasm import: ' + module + '.' + name);
const runtime = await createFormatsCpu(binary, {printErr() {}});
function boozerSpectrum(bytes) {
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  let offset = 12;
  for (let i = 0; i < 3; ++i) offset += 4 + view.getInt32(offset, true);
  const metadata = Array.from({length: 13}, (_, i) => view.getInt32(offset + 4 * i, true));
  offset += 13 * 4;
  const lasym = view.getInt32(offset, true); offset += 12;
  const surfaces = metadata[1] - metadata[7], modes = (metadata[10] + 1) * (metadata[11] + 1);
  offset += 16 * surfaces + 8 * modes;
  const begin = offset, count = surfaces * modes;
  const fields = Array.from({length: lasym ? 12 : 6}, () => Array.from({length: count}, () => {
    const value = view.getFloat64(offset, true); offset += 8; return value;
  }));
  return {fields, prefix: bytes.slice(0, begin), suffix: bytes.slice(offset)};
}
assert(runtime.cpu.memory.buffer instanceof ArrayBuffer, 'file conversion uses unshared memory');
assert.equal(fixture.exports.equilibrium_test_netcdf4(), 0, 'compressed HDF5-backed NetCDF-4 must round-trip');
assert.throws(() => runtime.convert(new Uint8Array(3), 'h5'), /truncated/);
assert.throws(() => runtime.convert(new Uint8Array(8), 'zip'), /Unknown/);
assert.equal(fixture.exports.equilibrium_test_format_fixture(6), 0);
assert.throws(() => runtime.convert(files.read('/fixture.bin'), 'h5'), /raw and folded input/);
files.remove('/fixture.bin');
if (output) await mkdir(output, {recursive: true});
for (const kind of [0, 1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 12]) {
  assert.equal(fixture.exports.equilibrium_test_format_fixture(kind), 0, 'create fixture ' + kind);
  const original = files.read('/fixture.bin'); files.remove('/fixture.bin');
  if (kind >= 4 && kind !== 9 && kind !== 10) {
    const before = boozerSpectrum(original), after = boozerSpectrum(runtime.convert(original, 'bin'));
    assert.deepEqual(after.prefix, before.prefix, 'Boozer metadata and profiles remain exact');
    assert.deepEqual(after.suffix, before.suffix, 'Boozer sampled fields, tracing and angle maps remain exact');
    // Each parity pair passes through scaled sums/differences in the existing
    // library. Bound rounding by both inputs, including cancellation cases.
    for (let field = 0; field < before.fields.length; field += 2)
      for (let i = 0; i < before.fields[field].length; ++i) {
        const bound = 8 * Number.EPSILON * (Math.abs(before.fields[field][i]) + Math.abs(before.fields[field + 1][i]));
        for (let j = field; j <= field + 1; ++j)
          assert(Math.abs(after.fields[j][i] - before.fields[j][i]) <= bound,
            'Boozer spectral re-encoding exceeded its parity-pair rounding bound');
      }
  }
  if (output) await writeFile(join(output, kind + '.bin'), original);
  for (const format of ['nc', 'h5']) {
    const converted = runtime.convert(original, format);
    assert(converted.length > 100);
    assert.equal(Buffer.from(converted.subarray(0, format === 'nc' ? 3 : 8)).toString('hex'),
      format === 'nc' ? '434446' : '894844460d0a1a0a');
    if (kind < 4 || kind === 9 || kind === 10) assert.deepEqual(runtime.convert(converted, 'bin', format), original,
      `fixture ${kind} preserves every binary field and provenance through ${format}`);
    if ((kind < 4 || kind === 9 || kind === 10) && format === 'nc') {
      files.write('/fixture.nc', converted);
      assert.equal(fixture.exports.equilibrium_test_netcdf4_copy(), 0);
      const compressed = files.read('/fixture4.nc');
      files.remove('/fixture.nc'); files.remove('/fixture4.nc');
      assert.deepEqual(runtime.convert(compressed, 'bin', 'nc'), original,
        'compressed NetCDF-4 retains the complete original state and report');
    }
    assert.throws(() => runtime.files.read('/source.bin'), /could not open/);
    assert.throws(() => runtime.files.read('/converted.' + format), /could not open/);
    if (output) await writeFile(join(output, kind + '.' + format), converted);
  }
}
console.log('PASS: standalone formats, compressed NetCDF-4, exact VMEC state/provenance round trips, asymmetric/free metadata, Boozer v6/v7/v8 exports, cleanup and recovered errors');
