// Central data-loading helpers. All YAML under ../data (generated) and
// ../config (curated) is read at build time.
import { readFileSync } from 'node:fs';
import { join } from 'node:path';
import YAML from 'yaml';

// At build/dev time process.cwd() is the astro project dir (website/astro).
const PROJECT_DIR = process.cwd();
const DATA_DIR = join(PROJECT_DIR, '..', 'data');     // website/data (generated)
const CONFIG_DIR = join(PROJECT_DIR, '..', 'config');  // website/config (curated)

export function loadYaml<T = any>(name: string): T {
  const raw = readFileSync(join(DATA_DIR, name), 'utf-8');
  return YAML.parse(raw) as T;
}

export function loadConfig<T = any>(name: string): T {
  const raw = readFileSync(join(CONFIG_DIR, name), 'utf-8');
  return YAML.parse(raw) as T;
}

export interface EdamTerm {
  uri: string;
  label: string;
}

export interface ToolIndexEntry {
  id: string;
  name: string;
  version: string;
  description: string;
  edam_operations: EdamTerm[];
  edam_topics: EdamTerm[];
  biotools: string | null;
  doi: string | null;
  panel: string | null;
  panel_id: string | null;
  owner: string;
  repo: string;
  path: string;
  tests: number;
  updated: string;
  input_types: string[];
  output_types: string[];
}

export interface ToolFull extends ToolIndexEntry {
  profile: string;
  tool_type: string;
  license: string | null;
  hidden: boolean;
  icon: string | null;
  xrefs: { value: string; type: string }[];
  panel_section_id: string | null;
  panel_section_name: string | null;
  source_path: string;
  tool_dir: string;
  test_count: number;
  requirements: any[];
  container_requirements: any[];
  container_links: {
    conda: { name: string; version: string; url: string }[];
    docker: string | null;
    singularity: string | null;
  };
  creators: any[];
  help: string;
  help_html: string;
  help_format: string | null;
  inputs: any[];
  outputs: any[];
  tags: string[];
  updated: string;
  source_dependency_hash?: string;
  extracted_at?: string;
}

export interface Contributor {
  id: string;
  name: string;
  github: string | null;
  orcid: string | null;
  email: string | null;
  url: string | null;
  tools: string[];
  is_member: boolean;
  tool_count: number;
  commits: number;
  first_commit: string | null;
  gtn_handle: string | null;
}

export interface GtnContributor {
  name: string;
  email?: string;
  orcid?: string;
  bio?: string;
  affiliations?: string[];
  former_affiliations?: string[];
  location?: { country: string; lat: number; lon: number };
  joined?: string;
  fediverse?: string;
  bluesky?: string;
  linkedin?: string;
  matrix?: string;
  twitter?: string;
  elixir_node?: string;
  url?: string;
}

export interface GtnOrganisation {
  name: string;
  short_name?: string;
  url?: string;
  avatar?: string;
  ror?: string;
  bio?: string;
}

interface GtnPeopleData {
  fetched_at: string;
  contributors: Record<string, GtnContributor>;
  organisations: Record<string, GtnOrganisation>;
}

export interface Organisation {
  id: string;
  name: string;
  url: string | null;
  identifier: string | null;
  tools: string[];
  tool_count: number;
}

export interface ServerConfig {
  servers: { name: string; url: string }[];
  toolshed: { base: string; host: string; view: string };
  iuc: { repo: string; issues: string; guides: string };
  ecosystem: {
    galaxy_hub: string;
    galaxy_hub_iuc: string;
    galaxy_tools: string;
    iwc: string;
    chat: string;
  };
}

export interface ToolAvailability {
  servers: Record<string, {
    url: string;
    tool_ids: string[];
    repos: string[];
    tool_count: number;
    repo_count: number;
  }>;
}

export const siteConfig = loadConfig<ServerConfig>('site.yaml');

// GTN hall-of-fame URL prefixes for cross-linking.
const GTN_HOF_URL = 'https://training.galaxyproject.org/training-material/hall-of-fame';
const HUB_HOF_URL = 'https://galaxyproject.org/hall-of-fame';

let _gtn: GtnPeopleData | null = null;

export function getGtnPeople(): GtnPeopleData {
  if (!_gtn) {
    try {
      _gtn = loadYaml<GtnPeopleData>('gtn_people.yaml');
    } catch {
      _gtn = { fetched_at: '', contributors: {}, organisations: {} };
    }
  }
  return _gtn;
}

/** Get GTN contributor metadata by handle, or null if not in GTN. */
export function getGtnContributor(handle: string | null): GtnContributor | null {
  if (!handle) return null;
  return getGtnPeople().contributors[handle] ?? null;
}

/** Get the GTN hall-of-fame URL for a handle. */
export function gtnHallOfFameUrl(handle: string | null): string | null {
  return handle ? `${GTN_HOF_URL}/${handle}/` : null;
}

/** Get the Galaxy Hub hall-of-fame URL for a handle. */
export function hubHallOfFameUrl(handle: string | null): string | null {
  return handle ? `${HUB_HOF_URL}/${handle}/` : null;
}

/** Get the URL for an affiliation by its handle from GTN organisations. */
export function affiliationUrl(affiliation: string): string | null {
  return getGtnPeople().organisations[affiliation]?.url ?? null;
}

let _availability: ToolAvailability | null = null;

export function getToolAvailability(): ToolAvailability {
  if (!_availability) {
    try {
      _availability = loadYaml<ToolAvailability>('tool_availability.yaml');
    } catch {
      _availability = { servers: {} };
    }
  }
  return _availability;
}

export function isToolInstalled(serverName: string, owner: string, repo: string, toolId: string): boolean {
  const avail = getToolAvailability();
  const srv = avail.servers[serverName];
  if (!srv) return false;
  // Match by tool short ID first (most reliable — handles cases where the
  // ToolShed repo name differs from the IUC .shed.yml repo name).
  if (srv.tool_ids?.includes(toolId)) return true;
  // Fall back to owner/repo matching.
  return srv.repos?.includes(`${owner}/${repo}`) ?? false;
}

export function toolShedUrl(owner: string, repo: string): string {
  return `${siteConfig.toolshed.base}${siteConfig.toolshed.view}${owner}/${repo}`;
}

export function runOnUrl(server: { url: string }, owner: string, repo: string, id: string, version: string): string {
  const shedHost = siteConfig.toolshed.host;
  const guid = `${shedHost}/repos/${owner}/${repo}/${id}/${version}`;
  return `${server.url}/?tool_id=${encodeURIComponent(guid)}`;
}

// Read a per-tool YAML file by owner/repo/id.
export function loadTool(owner: string, repo: string, id: string): ToolFull {
  const raw = readFileSync(join(DATA_DIR, 'tools', owner, repo, `${id}.yaml`), 'utf-8');
  return YAML.parse(raw) as ToolFull;
}
