// Compare the serialized Fourier curve with points sampled in the original
// Miller angle. This check is independent of the generator's Fourier angle.
import assert from 'node:assert/strict';
import { execFileSync } from 'node:child_process';
import { mkdtempSync, readFileSync, rmSync, writeFileSync } from 'node:fs';
import { join } from 'node:path';

const scratch = mkdtempSync(join(process.cwd(), 'miller-cli-'));
try {
  const input = join(scratch, 'source.g');
  const output = join(scratch, 'batch');
  const gfile = 'COCOS=5'.padEnd(48) + '   0   8   8\n' +
    '2 2 3 2 0.25\n3 0.25 0 1 2\n1000 0 0 0 0\n0 0 0 0 0\n' +
    [6, 400, 0, 0].map(x => Array(8).fill(x).join(' ')).join('\n') + '\n' +
    Array(64).fill(0).join(' ') + '\n' + Array(8).fill(-2).join(' ') + '\n' +
    '4 0\n3.6 0.25 3 1.09 2.4 0.25 3 -0.59\n';
  writeFileSync(input, gfile);
  execFileSync(process.argv[2], [input, '--miller-batch', output, '--mpol', '16',
    '--r0', '3.3', '--kappa-upper', '1.6', '--kappa-lower', '1.2',
    '--delta-upper', '0.2', '--delta-lower', '-0.2']);
  const p = JSON.parse(readFileSync(join(output, 'case_0.json'), 'utf8'));
  assert.equal(p.lasym, true);
  assert.deepEqual(p.raxis_c, [3.3]);
  assert.deepEqual(p.zaxis_c, [0.25]);
  const sum = (modes, fn, t) =>
    modes.reduce((v, {m, value}) => v + value * fn(m*t), 0);
  const polygon = Array.from({length: 4096}, (_, i) => {
    const t = 2*Math.PI*i/4096;
    return [sum(p.rbc, Math.cos, t) + sum(p.rbs, Math.sin, t),
            sum(p.zbs, Math.sin, t) + sum(p.zbc, Math.cos, t)];
  });
  for (let i = 0; i < 193; ++i) {
    const u = 2*Math.PI*(i + 0.37)/193;
    const upper = u < Math.PI;
    const r = 3.3 + 0.6*Math.cos(u + Math.asin(upper ? 0.2 : -0.2)*Math.sin(u));
    const z = 0.25 + 0.6*(upper ? 1.6 : 1.2)*Math.sin(u);
    let best = Infinity;
    for (let j = 0; j < polygon.length; ++j) {
      const [x,y] = polygon[j], [x1,y1] = polygon[(j+1)%polygon.length];
      const dx = x1-x, dy = y1-y;
      const t = Math.max(0, Math.min(1, ((r-x)*dx+(z-y)*dy)/(dx*dx+dy*dy)));
      best = Math.min(best, Math.hypot(r-x-t*dx, z-y-t*dy));
    }
    assert.ok(best <= 0.6e-3, `Miller arc distance ${best} m at u=${u}`);
  }
} finally {
  rmSync(scratch, {recursive: true, force: true});
}
