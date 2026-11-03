import assert from 'node:assert/strict';
import test from 'node:test';

import { searchStateFromUrl, toolMatchesSearch } from './toolSearch.mjs';

const tool = {
  name: 'bwa mem',
  id: 'bwa_mem',
  description: 'Burrows-Wheeler aligner for short reads',
  inputTypes: ['fastq', 'fasta'],
  outputTypes: ['bam'],
};

test('searchStateFromUrl reads search and filter URLs', () => {
  assert.deepEqual(searchStateFromUrl('?q=bwa&input=fastq&output=bam'), {
    query: 'bwa',
    inputType: 'fastq',
    outputType: 'bam',
  });
});

test('toolMatchesSearch matches name, identifier, and description queries', () => {
  assert.equal(toolMatchesSearch(tool, { query: 'BWA' }), true);
  assert.equal(toolMatchesSearch(tool, { query: 'bwa_mem' }), true);
  assert.equal(toolMatchesSearch(tool, { query: 'short reads' }), true);
  assert.equal(toolMatchesSearch(tool, { query: 'variant caller' }), false);
});

test('toolMatchesSearch combines query with datatype filters', () => {
  assert.equal(toolMatchesSearch(tool, { query: 'aligner', inputType: 'fastq', outputType: 'bam' }), true);
  assert.equal(toolMatchesSearch(tool, { query: 'aligner', inputType: 'txt', outputType: 'bam' }), false);
  assert.equal(toolMatchesSearch(tool, { query: 'aligner', inputType: 'fastq', outputType: 'vcf' }), false);
});
