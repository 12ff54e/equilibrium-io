// cuMES result binaries v10/v11. The digest excludes headers and provenance.
function byteView(input) {
  if (ArrayBuffer.isView(input))
    return new Uint8Array(input.buffer, input.byteOffset, input.byteLength);
  if (input instanceof ArrayBuffer) return new Uint8Array(input);
  throw new TypeError('Expected binary equilibrium bytes');
}

function product(...values) {
  const size = values.reduce((a, b) => a * b, 1);
  if (!Number.isSafeInteger(size) || size < 0)
    throw new RangeError('Scientific payload size overflow');
  return size;
}

export function scientificPayloadRanges(input) {
  const bytes = byteView(input);
  if (bytes.length < 20) throw new RangeError('Truncated result header');
  const magic = 'CUMES001';
  if (![...magic].every((value, i) => bytes[i] === value.charCodeAt(0)))
    throw new Error('Not a cuMES result binary');
  const view = new DataView(bytes.buffer, bytes.byteOffset, bytes.byteLength);
  const version = view.getInt32(8, true);
  const ns = view.getInt32(12, true), mnmax = view.getInt32(16, true);
  if (![10, 11].includes(version) || ns < 2 || mnmax < 1)
    throw new Error(`Unsupported output shape/version: ${version}/${ns}/${mnmax}`);
  const start = 24;
  if (bytes.length < start) throw new RangeError('Truncated component count');
  const components = view.getInt32(20, true);
  if (![6, 12].includes(components)) throw new Error('Invalid spectral component count');
  const spectralEnd = start + product(components, ns, mnmax, 8);
  if (!Number.isSafeInteger(spectralEnd + 8) || spectralEnd + 8 > bytes.length)
    throw new RangeError('Truncated spectral payload');
  const ntheta = view.getInt32(spectralEnd, true);
  const nzeta = view.getInt32(spectralEnd + 4, true);
  if (ntheta < 1 || nzeta < 1) throw new Error('Missing derived output fields');
  const end = spectralEnd + 8 + product(7 * (ns - 1) + 6 * ns, ntheta, nzeta, 8);
  if (!Number.isSafeInteger(end) || end > bytes.length)
    throw new RangeError('Truncated scientific output');
  return {version, ns, mnmax, components, ntheta, nzeta, scientificBytes: end,
    stateRange: [start, spectralEnd], fieldRange: [spectralEnd + 8, end]};
}

export async function digestScientificPayload(input) {
  const bytes = byteView(input);
  const {stateRange, fieldRange, components, ...metadata} = scientificPayloadRanges(bytes);
  const digest = async ([begin, end]) => Array.from(new Uint8Array(
    await globalThis.crypto.subtle.digest('SHA-256', bytes.subarray(begin, end))))
    .map(value => value.toString(16).padStart(2, '0')).join('');
  return {...metadata, state: await digest(stateRange), fields: await digest(fieldRange)};
}
