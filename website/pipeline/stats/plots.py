from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import seaborn as sns
import pandas as pd

from .data_model import ToolEvent

sns.set_style('whitegrid')
plt.rcParams['figure.dpi'] = 150
plt.rcParams['font.size'] = 10


def _monthly_counts(
    events: list[ToolEvent],
) -> pd.DataFrame:
    counts: dict[tuple[int, int], int] = defaultdict(int)
    for e in events:
        if e.commit_date is not None:
            key = (e.commit_date.year, e.commit_date.month)
            counts[key] += 1

    df = pd.DataFrame([
        {'year': y, 'month': m, 'count': c}
        for (y, m), c in sorted(counts.items())
    ])
    if df.empty:
        return df
    df['date'] = pd.to_datetime(
        df['year'].astype(str) + '-' + df['month'].astype(str) + '-01'
    )
    return df.sort_values('date')


def _cumulative_tools(events: list[ToolEvent]) -> pd.DataFrame:
    new_events = sorted(
        [e for e in events if e.event_type == 'new'],
        key=lambda e: e.commit_date or datetime.min,
    )
    rows = []
    cumulative = 0
    for e in new_events:
        cumulative += 1
        rows.append({
            'date': e.commit_date,
            'count': cumulative,
        })
    return pd.DataFrame(rows)


def plot_new_tools_per_month(
    events: list[ToolEvent],
    output_dir: str,
):
    new_tools = [e for e in events if e.event_type == 'new']
    df = _monthly_counts(new_tools)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(df['date'], df['count'], width=20, color='#2e86ab', edgecolor='white')
    ax.set_title('New Galaxy Tools Added Per Month', fontsize=14, fontweight='bold')
    ax.set_xlabel('Date')
    ax.set_ylabel('New Tools')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'new_tools_per_month.png')
    plt.close(fig)


def plot_updates_per_month(
    events: list[ToolEvent],
    output_dir: str,
):
    updates = [e for e in events if e.event_type == 'update']
    df = _monthly_counts(updates)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.bar(df['date'], df['count'], width=20, color='#ca6702', edgecolor='white')
    ax.set_title('Tool Version Updates Per Month', fontsize=14, fontweight='bold')
    ax.set_xlabel('Date')
    ax.set_ylabel('Updates')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'updates_per_month.png')
    plt.close(fig)


def plot_combined_monthly(
    events: list[ToolEvent],
    output_dir: str,
):
    new_df = _monthly_counts([e for e in events if e.event_type == 'new'])
    up_df = _monthly_counts([e for e in events if e.event_type == 'update'])

    merged = pd.merge(
        new_df[['date', 'count']].rename(columns={'count': 'new'}),
        up_df[['date', 'count']].rename(columns={'count': 'update'}),
        on='date', how='outer',
    ).fillna(0).sort_values('date')

    if merged.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 6))
    x = merged['date']
    width = 20
    ax.bar(x, merged['new'], width, label='New Tools', color='#2e86ab')
    ax.bar(x, merged['update'], width, label='Updates',
           bottom=merged['new'], color='#ca6702')
    ax.set_title('Monthly Galaxy Tool Activity', fontsize=14, fontweight='bold')
    ax.set_xlabel('Date')
    ax.set_ylabel('Count')
    ax.legend()
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'combined_monthly.png')
    plt.close(fig)


def plot_cumulative_tools(
    events: list[ToolEvent],
    output_dir: str,
):
    df = _cumulative_tools(events)
    if df.empty:
        return

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.fill_between(df['date'], df['count'], alpha=0.3, color='#2e86ab')
    ax.plot(df['date'], df['count'], color='#2e86ab', linewidth=1.5)
    ax.set_title('Cumulative Galaxy Tools in Repository', fontsize=14, fontweight='bold')
    ax.set_xlabel('Date')
    ax.set_ylabel('Total Tools')
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'cumulative_tools.png')
    plt.close(fig)


def plot_top_updated_tools(
    events: list[ToolEvent],
    output_dir: str,
    top_n: int = 20,
):
    update_counts: Counter = Counter()
    for e in events:
        if e.event_type == 'update':
            update_counts[e.tool_id] += 1

    if not update_counts:
        return

    top = update_counts.most_common(top_n)
    names, counts = zip(*reversed(top))

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.barh(names, counts, color='#2e86ab', edgecolor='white')
    ax.set_title(
        f'Top {top_n} Most Updated Tools', fontsize=14, fontweight='bold'
    )
    ax.set_xlabel('Number of Updates')
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'top_updated_tools.png')
    plt.close(fig)


def plot_activity_heatmap(
    events: list[ToolEvent],
    output_dir: str,
):
    heat: dict[tuple[int, int], int] = defaultdict(int)
    for e in events:
        if e.commit_date is not None:
            heat[(e.commit_date.year, e.commit_date.month)] += 1

    if not heat:
        return

    years = sorted({y for y, m in heat})
    months = list(range(1, 13))
    data = []
    for y in years:
        row = []
        for m in months:
            row.append(heat.get((y, m), 0))
        data.append(row)

    fig, ax = plt.subplots(figsize=(14, max(4, len(years) * 0.4)))
    sns.heatmap(
        data, ax=ax,
        xticklabels=['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                     'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'],
        yticklabels=years,
        cmap='YlOrRd',
        linewidths=0.5,
        cbar_kws={'label': 'Deployment Events'},
    )
    ax.set_title('Tool Activity Heatmap (Year × Month)', fontsize=14, fontweight='bold')
    ax.set_xlabel('Month')
    ax.set_ylabel('Year')
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'activity_heatmap.png')
    plt.close(fig)


def plot_change_type_breakdown(
    events: list[ToolEvent],
    output_dir: str,
):
    update_types: Counter = Counter()
    for e in events:
        if e.event_type == 'update':
            update_types[e.change_type or 'unknown'] += 1

    if not update_types:
        return

    labels_map = {
        'upstream': 'Upstream Version Bump',
        'wrapper': 'Wrapper Revision Bump',
        'version': 'Version Change',
        'initial': 'Initial Addition',
        None: 'Unknown',
        'unknown': 'Unknown',
    }

    labels = []
    sizes = []
    colors = ['#2e86ab', '#ca6702', '#52b788', '#e9c46a', '#adb5bd']
    for ct, cnt in update_types.most_common():
        labels.append(labels_map.get(ct, ct))
        sizes.append(cnt)

    fig, ax = plt.subplots(figsize=(8, 8))
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, autopct='%1.1f%%',
        colors=colors[:len(sizes)],
        startangle=90,
    )
    ax.set_title('Update Type Breakdown', fontsize=14, fontweight='bold')
    fig.tight_layout()
    fig.savefig(Path(output_dir) / 'change_type_breakdown.png')
    plt.close(fig)


def generate_all_plots(events: list[ToolEvent], output_dir: str):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    plot_new_tools_per_month(events, output_dir)
    plot_updates_per_month(events, output_dir)
    plot_combined_monthly(events, output_dir)
    plot_cumulative_tools(events, output_dir)
    plot_top_updated_tools(events, output_dir)
    plot_activity_heatmap(events, output_dir)
    plot_change_type_breakdown(events, output_dir)
