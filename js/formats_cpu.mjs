export async function createFormatsCpu(source, options = {}) {
  const load = path => {
    const url = new URL(path, import.meta.url);
    url.search = new URL(import.meta.url).search;
    return import(url.href);
  };
  const [{createCpuWasm}, {cpuFiles}] = await Promise.all([load('./cpu_wasm.mjs'), load('./cpu_files.mjs')]);
  const cpu = await createCpuWasm(source, options), files = cpuFiles(cpu);
  return {
    cpu, files,
    convert(bytes, format, inputFormat = 'bin') {
      if (!['bin', 'nc', 'h5'].includes(format) || !['bin', 'nc', 'h5'].includes(inputFormat))
        throw Error('Unknown scientific file format.');
      const input = '/source.' + inputFormat, output = '/converted.' + format;
      try {
        files.write(input, bytes);
        const code = cpu.withString(input, from => cpu.withString(output, to => cpu.exports.equilibrium_format_convert(from, to)));
        if (code) throw Error(cpu.string(cpu.exports.equilibrium_format_error()));
        return files.read(output);
      } finally {
        files.remove(input); files.remove(output);
      }
    }
  };
}
