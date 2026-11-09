import assert from 'node:assert/strict';
import {createHash, webcrypto} from 'node:crypto';
import {scientificPayloadRanges, digestScientificPayload} from '../js/binary_payload.mjs';

globalThis.crypto ??= webcrypto;

function fixture(version = 11, components = 6) {
  const start = version >= 9 ? 24 : 20;
  const stateSize = components * 2 * 1 * 8, fieldSize = 19 * 3 * 2 * 8;
  const bytes = new Uint8Array(start + stateSize + 8 + fieldSize + 16);
  bytes.set(new TextEncoder().encode('CUMES001'));
  const view = new DataView(bytes.buffer);
  view.setInt32(8, version, true);
  view.setInt32(12, 2, true);
  view.setInt32(16, 1, true);
  if (version >= 9) view.setInt32(20, components, true);
  bytes.fill(0xa5, start, start + stateSize);
  view.setInt32(start + stateSize, 3, true);
  view.setInt32(start + stateSize + 4, 2, true);
  bytes.fill(0x5a, start + stateSize + 8, bytes.length - 16);
  return {bytes, start, stateSize, fieldSize};
}
const hash = bytes => createHash('sha256').update(bytes).digest('hex');
for (const version of [10, 11]) {
  const {bytes, start, stateSize, fieldSize} = fixture(version, version >= 9 ? 12 : 6);
  const ranges = scientificPayloadRanges(bytes);
  assert.equal(ranges.components, version >= 9 ? 12 : 6);
  assert.deepEqual(ranges.stateRange, [start, start + stateSize]);
  assert.deepEqual(ranges.fieldRange, [start + stateSize + 8, bytes.length - 16]);
  const digest = await digestScientificPayload(bytes);
  assert(!Object.hasOwn(digest, 'components')); // Preserve the existing digest record.
  assert.equal(digest.state, hash(new Uint8Array(stateSize).fill(0xa5)));
  assert.equal(digest.fields, hash(new Uint8Array(fieldSize).fill(0x5a)));
  bytes.fill(0xff, bytes.length - 16); // Provenance cannot affect the digest.
  assert.deepEqual(await digestScientificPayload(bytes.buffer), digest);
  const padded = new Uint8Array(bytes.length + 11);
  padded.set(bytes, 7);
  assert.deepEqual(await digestScientificPayload(padded.subarray(7, 7 + bytes.length)), digest);
  for (const length of [0, 19, start - 1, start + stateSize - 1, bytes.length - 17])
    assert.throws(() => scientificPayloadRanges(bytes.subarray(0, length)));
}
for (const [offset, value] of [[8, 0], [8, 8], [8, 9], [8, 12], [12, 0], [16, -1], [20, 7]]) {
  const {bytes} = fixture(10, 12);
  new DataView(bytes.buffer).setInt32(offset, value, true);
  assert.throws(() => scientificPayloadRanges(bytes));
}
{
  const {bytes} = fixture();
  bytes[0] = 0;
  assert.throws(() => scientificPayloadRanges(bytes), /Not a cuMES/);
}
{
  const {bytes} = fixture();
  const view = new DataView(bytes.buffer);
  view.setInt32(12, 0x7fffffff, true);
  view.setInt32(16, 0x7fffffff, true);
  assert.throws(() => scientificPayloadRanges(bytes), /overflow/);
}
{
  const {bytes, start, stateSize} = fixture();
  const view = new DataView(bytes.buffer);
  view.setInt32(start + stateSize, 0, true);
  assert.throws(() => scientificPayloadRanges(bytes), /Missing derived/);
  view.setInt32(start + stateSize, 0x7fffffff, true);
  view.setInt32(start + stateSize + 4, 0x7fffffff, true);
  assert.throws(() => scientificPayloadRanges(bytes), /overflow/);
}
console.log('binary payload ranges and digests pass');
