"""
Galaxy Toolshed Deployment Analyzer for tools-iuc.

Walks the entire git history of the tools-iuc repository along the
origin/main first-parent chain, discovers every <tool id="..."> across
all XML files in tools/, resolves version macros (TOOL_VERSION,
VERSION_SUFFIX, WRAPPER_VERSION, etc.) at each commit, and classifies
every version change as:

  - new       — first appearance of a tool id on origin/main
  - upstream  — the upstream tool version (TOOL_VERSION / VERSION) changed
  - wrapper   — the Galaxy wrapper revision (VERSION_SUFFIX / WRAPPER_VERSION) changed
  - version   — version changed via custom tokens not in our known lists
                (e.g. @STACKS_VERSION@, @galaxy_version@) or hardcoded strings

Outputs raw data (CSV + JSON), a summary table, and 7 matplotlib/seaborn
plots to --output-dir (default: .stats/).

Usage:
    python -m stats.main --repo /path/to/tools-iuc --verbose
"""

import argparse
import sys
import time
from pathlib import Path

from .discovery import discover_tools
from .git_analyzer import analyze_tools
from .output import write_events_csv, write_events_json, write_summary
from .plots import generate_all_plots


def main():
    parser = argparse.ArgumentParser(
        description='Analyze Galaxy tool additions and updates from git history.'
    )
    parser.add_argument(
        '--repo', default='.',
        help='Path to the tools-iuc git repository (default: current dir)',
    )
    parser.add_argument(
        '--branch', default=None,
        help='Git branch to analyze (default: auto-detect origin/main or master)',
    )
    parser.add_argument(
        '--tools-path', default='tools',
        help='Path to tools subdirectory within repo (default: tools)',
    )
    parser.add_argument(
        '--output-dir', default='./.stats',
        help='Directory for output files (default: ./.stats)',
    )
    parser.add_argument(
        '--skip-plots', action='store_true',
        help='Skip plot generation',
    )
    parser.add_argument(
        '--format', choices=['csv', 'json', 'both'], default='both',
        help='Output format for raw data (default: both)',
    )
    parser.add_argument(
        '--verbose', '-v', action='store_true',
        help='Verbose output',
    )
    parser.add_argument(
        '--workers', type=int, default=1,
        help='Parallel workers for tool analysis (default: 1)',
    )

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()

    if args.verbose:
        print('Discovering tools...', file=sys.stderr)

    tools = discover_tools(repo_path=args.repo, tools_path=args.tools_path)
    if args.verbose:
        print(f'Found {len(tools)} Galaxy tools with <tool id=> in XML files',
              file=sys.stderr)

    if not tools:
        print('No tools found. Check --repo and --tools-path.', file=sys.stderr)
        sys.exit(1)

    t1 = time.time()
    if args.verbose:
        print(f'Discovery took {t1 - t0:.1f}s', file=sys.stderr)
        print('Analyzing git history...', file=sys.stderr)

    events = analyze_tools(
        repo_path=args.repo,
        tools=tools,
        branch=args.branch,
        tools_path=args.tools_path,
        verbose=args.verbose,
        workers=args.workers,
    )

    t2 = time.time()
    if args.verbose:
        print(f'Analysis took {t2 - t1:.1f}s', file=sys.stderr)
        print(f'Total events: {len(events)}', file=sys.stderr)

    fmt = args.format
    if fmt in ('csv', 'both'):
        csv_path = output_dir / 'tool_events.csv'
        write_events_csv(events, str(csv_path))
        if args.verbose:
            print(f'Wrote {csv_path}', file=sys.stderr)

    if fmt in ('json', 'both'):
        json_path = output_dir / 'tool_events.json'
        write_events_json(events, str(json_path))
        if args.verbose:
            print(f'Wrote {json_path}', file=sys.stderr)

    summary_path = output_dir / 'summary.txt'
    write_summary(events, str(summary_path))
    if args.verbose:
        print(f'Wrote {summary_path}', file=sys.stderr)

    if not args.skip_plots:
        if args.verbose:
            print('Generating plots...', file=sys.stderr)
        generate_all_plots(events, str(output_dir))
        if args.verbose:
            print(f'Plots saved to {output_dir}/', file=sys.stderr)

    t3 = time.time()
    if args.verbose:
        print(f'Total time: {t3 - t0:.1f}s', file=sys.stderr)

    additions = len([e for e in events if e.event_type == 'new'])
    updates = len([e for e in events if e.event_type == 'update'])
    print(f'Done. {additions} new tools, {updates} updates found.')
    print(f'Output in {output_dir}/')


if __name__ == '__main__':
    main()
