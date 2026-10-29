export function searchStateFromUrl(search) {
  const params = new URLSearchParams(search);
  return {
    query: params.get('q') ?? '',
    inputType: params.get('input') ?? '',
    outputType: params.get('output') ?? '',
  };
}

export function toolMatchesSearch(tool, state) {
  const query = (state.query || '').toLowerCase();
  const inputType = state.inputType || '';
  const outputType = state.outputType || '';
  const inputs = tool.inputTypes || [];
  const outputs = tool.outputTypes || [];
  const searchableText = (tool.searchableText || `${tool.name || ''} ${tool.id || ''} ${tool.description || ''}`).toLowerCase();

  return (!query || searchableText.includes(query))
    && (!inputType || inputs.includes(inputType))
    && (!outputType || outputs.includes(outputType));
}
