import { loadYaml } from '../lib/data';

type ToolShedRepo = {
  name?: string;
  owner?: string;
  description?: string;
};

function normalize(value: string | undefined): string {
  return (value || '').toLowerCase().replace(/[^a-z0-9]+/g, '_').replace(/^_+|_+$/g, '');
}

export function GET() {
  const index = loadYaml<{ repos: ToolShedRepo[] }>('toolshed_index.yaml');
  const repos = (index.repos || [])
    .filter((repo) => repo.name && repo.owner)
    .map((repo) => [
      normalize(repo.name),
      repo.owner,
      repo.name,
      (repo.description || '').slice(0, 120),
    ]);

  return new Response(JSON.stringify(repos), {
    headers: {
      'Content-Type': 'application/json; charset=utf-8',
      'Cache-Control': 'public, max-age=3600',
    },
  });
}
