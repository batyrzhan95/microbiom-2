"""Reproduce Figures 2 and 3 from the study extraction Excel workbook.

Google Colab: paste this script into a code cell and run it; select your .xlsx.
Local Python: python semen_microbiome_figures.py 'Study_Data_Extraction(1).xlsx'
Requires pandas, numpy, matplotlib and openpyxl. No meta-analysis is fitted.
"""

from pathlib import Path
import hashlib
import json
import platform
import sys
import zipfile

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import openpyxl

TEAL, MAGENTA = '#247D7B', '#A34770'
STYLE = {
    'font.family': 'DejaVu Sans', 'font.size': 10.5,
    'axes.labelcolor': '#263740', 'text.color': '#263740',
    'xtick.color': '#263740', 'ytick.color': '#263740',
    'pdf.fonttype': 42, 'svg.fonttype': 'none',
    'savefig.facecolor': 'white', 'figure.facecolor': 'white',
}

CAPTION_2 = (
    'Figure 2. Reported denominators for microbiome and sperm DNA analyses. '
    'Bars show study-specific microbiome sample sizes and the available sperm DNA '
    'fragmentation (SDF) summary denominator; these do not establish the number of '
    'paired microbiome–SDF observations. Alkan et al. (2025): 16 microbiome and '
    '14 TUNEL specimens, with correlation denominators varying by taxon. '
    'Garcia-Segura et al. (2022): 56 participants; the separate assay-specific SDF '
    'denominator is not reported (NR). He et al. (2024): the main 87-person '
    'microbiome cohort is shown; total enrollment was reported as 102. '
    'Mowla et al. (2025, Version of Record): the DNA-fragmentation categories in '
    'Table 1 total 196, whereas Appendix 3 assay counts total 223. This discrepancy '
    'does not establish that 27 participants were untested. Denominators are not pooled.'
)
CAPTION_3 = (
    'Figure 3. Source-specific correlations between seminal microbiome measures '
    'and sperm DNA outcomes. Panel A displays all 12 DNA fragmentation index '
    '(DFI) coefficients in Alkan et al. (2025), Table 5. Panel B displays all five '
    'DNA-damage outcomes in the Shannon diversity block of Garcia-Segura et al. '
    '(2022), Table 3, including Olive tail moment (OTM). Points are reported '
    'Spearman correlations. Horizontal stems connect zero to each coefficient '
    'and are not confidence intervals. Accompanying n values are pairwise '
    'denominators where explicitly reported; NR denotes not reported. Alkan '
    'Table 5 is itself a selected published table, rather than a complete inventory '
    'of tests. The panels use different exposures and assays and are not pooled. '
    'No confidence intervals or significance colors are added. Nominal '
    'associations do not establish multiplicity-adjusted evidence or causation.'
)


def numeric_column(frame, column, sheet, integer=False, lower=None, upper=None):
    """Reject invalid values without turning missing observations into zero."""
    original = frame[column]
    missing = original.isna() | original.astype(str).str.strip().str.lower().isin(
        ['', 'na', 'n/a', 'nr', 'not reported', 'none', 'nan']
    )
    values = pd.to_numeric(original.mask(missing), errors='coerce')
    bad = (~missing & values.isna()) | (values.notna() & ~np.isfinite(values))
    if integer:
        bad |= values.notna() & (values % 1 != 0)
    if lower is not None:
        bad |= values.notna() & (values < lower)
    if upper is not None:
        bad |= values.notna() & (values > upper)
    if bad.any():
        rows = frame.loc[bad, 'excel_row'].astype(int).tolist()
        raise ValueError(f'{sheet}.{column}: invalid value(s) at Excel row(s) {rows}.')
    frame[column] = values.astype(float)


def read_workbook(workbook):
    """Read the two plotting sheets; preserve workbook row order and provenance."""
    required = {
        'Samples': ['study_id', 'study_label', 'total_reported_n', 'microbiome_n',
                    'sdf_n', 'sdf_n_status', 'doi', 'source'],
        'Correlations': ['study_id', 'study_label', 'exposure', 'outcome', 'rho',
                         'pairwise_n', 'source_table', 'source_page', 'doi'],
    }
    data = {}
    with pd.ExcelFile(workbook, engine='openpyxl') as book:
        for sheet, columns in required.items():
            if sheet not in book.sheet_names:
                raise ValueError(f'Missing worksheet: {sheet}. Use the study extraction workbook.')
            frame = pd.read_excel(book, sheet_name=sheet)
            frame.columns = frame.columns.astype(str).str.strip()
            absent = set(columns) - set(frame.columns)
            if absent:
                raise ValueError(f'{sheet}: missing columns {sorted(absent)}. Headers must be in row 1.')
            frame['excel_row'] = np.arange(len(frame)) + 2
            frame = frame.dropna(subset=columns, how='all').copy()
            for col in ['study_id', 'study_label']:
                if frame[col].isna().any():
                    raise ValueError(f'{sheet}: {col} contains a blank identifier.')
                frame[col] = frame[col].astype(str).str.strip()
            data[sheet] = frame
    studies, correlations = data['Samples'], data['Correlations']
    expected_ids = {'Alkan_2025', 'Garcia_Segura_2022', 'He_2024', 'Mowla_2025'}
    if set(studies.study_id) != expected_ids or studies.study_id.duplicated().any():
        raise ValueError('Samples must contain one row for each of the four studies in this review.')
    for col in ['total_reported_n', 'microbiome_n', 'sdf_n']:
        numeric_column(studies, col, 'Samples', integer=True, lower=0)
    for col in ['microbiome_n', 'sdf_n']:
        if (studies[col] > studies.total_reported_n).any():
            raise ValueError(f'Samples: {col} exceeds total_reported_n.')
    numeric_column(correlations, 'rho', 'Correlations', lower=-1, upper=1)
    numeric_column(correlations, 'pairwise_n', 'Correlations', integer=True, lower=2)
    source = correlations.source_table.fillna('').astype(str).str.strip()
    outcome = correlations.outcome.fillna('').astype(str)
    exposure = correlations.exposure.fillna('').astype(str)
    alkan = correlations.loc[
        correlations.study_id.eq('Alkan_2025') & source.str.fullmatch(r'Table\s+5', case=False)
        & outcome.str.contains('TUNEL|DFI', case=False)
    ].copy()
    garcia = correlations.loc[
        correlations.study_id.eq('Garcia_Segura_2022') & source.str.fullmatch(r'Table\s+3', case=False)
        & exposure.str.contains('Shannon', case=False)
        & outcome.str.contains('TUNEL|Comet|OTM', case=False)
    ].copy()
    for label, panel, expected in [('Alkan Table 5', alkan, 12), ('Garcia-Segura Table 3', garcia, 5)]:
        if len(panel) != expected:
            raise ValueError(f'{label}: expected {expected} plotted rows, found {len(panel)}. Check extraction completeness.')
        if panel.duplicated(['exposure', 'outcome']).any():
            raise ValueError(f'{label}: duplicate exposure/outcome rows.')
    return studies, alkan, garcia


def draw_coverage(studies):
    """Figure 2: microbiome and SDF summary denominators, including explicit NR."""
    with plt.rc_context(STYLE):
        fig, ax = plt.subplots(figsize=(11.8, 7.6))
        fig.subplots_adjust(left=.265, right=.94, top=.76, bottom=.29)
        y = np.arange(len(studies))
        for column, offset, color, label in [
            ('microbiome_n', -.16, TEAL, 'Microbiome analysis'),
            ('sdf_n', .16, MAGENTA, 'SDF summary denominator'),
        ]:
            values = studies[column].to_numpy(dtype=float)
            present = np.isfinite(values)
            ax.barh(y[present] + offset, values[present], height=.27, color=color, label=label)
            for i, value in enumerate(values):
                ax.text(4 if np.isnan(value) else value + 4, y[i] + offset,
                        'NR' if np.isnan(value) else str(int(value)), va='center')
        ax.set_yticks(y, studies.study_label.tolist())
        ax.set_ylim(len(studies) - .5, -.5)
        largest = np.nanmax(studies[['microbiome_n', 'sdf_n']].to_numpy(dtype=float))
        ax.set_xlim(0, max(largest * 1.12, 25))
        ax.set_xlabel('Number of participants / specimens reported', labelpad=14)
        ax.set_axisbelow(True)
        ax.grid(axis='x', color='#E5EBED', linewidth=.7)
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        ax.spines['bottom'].set_color('#819198')
        ax.tick_params(axis='y', length=0, pad=13)
        fig.suptitle('Figure 2. Reported analytic sample coverage', x=.04, y=.96,
                     ha='left', fontsize=16, weight='bold')
        fig.text(.04, .90, 'Study-specific microbiome and SDF summary denominators', fontsize=11.5)
        fig.legend(loc='upper left', bbox_to_anchor=(.257, .865), ncol=2, frameon=False, fontsize=10.5)
        fig.text(.04, .12,
                 'NR = not reported. Garcia-Segura: assay-specific SDF denominator not stated.\n'
                 'Alkan: pairwise n varies by taxon. He: main 87-person microbiome cohort only.',
                 linespacing=1.7, fontsize=10)
        fig.text(.04, .04,
                 'Mowla: Table 1 categorizes 196 men; Appendix 3 reports 223 SDF assays.\n'
                 'The difference is an unresolved reporting discrepancy; missing testing is not established.',
                 linespacing=1.6, fontsize=10)
        return fig


def draw_correlations(alkan, garcia):
    """Figure 3: all 12 Alkan DFI and 5 Garcia-Segura Shannon coefficients."""
    labels = {
        'TUNEL_percent': 'TUNEL-positive sperm (%)',
        'Alkaline_Comet_percent': 'Alkaline Comet (%)',
        'Alkaline_Comet_OTM': 'Alkaline Comet (OTM)',
        'Neutral_Comet_percent': 'Neutral Comet (%)',
        'Neutral_Comet_OTM': 'Neutral Comet (OTM)',
    }
    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(2, 1, figsize=(12.8, 12.5), gridspec_kw={'height_ratios': [12, 5]})
        fig.subplots_adjust(left=.25, right=.77, top=.845, bottom=.22, hspace=.64)
        for ax, data, color, heading, subtitle, panel_a in [
            (axes[0], alkan, TEAL, 'A. Alkan et al. (2025)',
             'Bacterial genera and TUNEL DFI; all 12 rows in Table 5', True),
            (axes[1], garcia, MAGENTA, 'B. Garcia-Segura et al. (2022)',
             'Shannon diversity and all five sperm DNA outcomes in Table 3', False),
        ]:
            y = np.arange(len(data))
            rho = data.rho.to_numpy(dtype=float)
            present = np.isfinite(rho)
            ax.hlines(y[present], 0, rho[present], color='#C7D1D5', linewidth=1.5)
            ax.scatter(rho[present], y[present], color=color, s=52, zorder=4)
            tick_labels = data.exposure.tolist() if panel_a else [labels.get(v, v.replace('_', ' ')) for v in data.outcome]
            ax.set_yticks(y, tick_labels)
            if panel_a:
                for tick in ax.get_yticklabels():
                    tick.set_fontstyle('italic')
            ax.set_ylim(len(data) - .35, -.65)
            ax.set_xlim(-1.06, 1.06)
            ax.set_xticks([-1, -.5, 0, .5, 1])
            ax.axvline(0, color='#6B7C83', linestyle='--', linewidth=.9)
            ax.set_axisbelow(True)
            ax.grid(axis='x', color='#E5EBED', linestyle=':')
            for spine in ['top', 'right', 'left']:
                ax.spines[spine].set_visible(False)
            ax.spines['bottom'].set_color('#819198')
            ax.tick_params(axis='y', length=0, pad=13)
            ax.set_xlabel('Reported Spearman correlation', labelpad=10)
            ax.set_title(heading, loc='left', fontweight='bold', fontsize=12, pad=39)
            ax.text(0, 1.065, subtitle, transform=ax.transAxes, fontsize=10.5)
            for yy, row in zip(y, data.itertuples(index=False)):
                n_text = 'NR' if pd.isna(row.pairwise_n) else str(int(row.pairwise_n))
                rho_text = 'NR' if pd.isna(row.rho) else f'{row.rho:+.3f}'
                ax.text(1.045, yy, f'rho = {rho_text}; n = {n_text}',
                        transform=ax.get_yaxis_transform(), va='center', fontsize=10)
        fig.suptitle('Figure 3. Reported microbiome–DNA correlations', x=.04, y=.974,
                     ha='left', fontsize=16, weight='bold')
        fig.text(.04, .933, 'Descriptive source-specific estimates; no pooled effect or confidence intervals', fontsize=11.5)
        fig.text(.04, .075,
                 'NR = not reported. n is the explicitly reported pairwise denominator.\n'
                 'Alkan Table 5 is a selected published table and does not include every tested genus.\n'
                 'OTM = Olive tail moment. Stems connect zero to each coefficient; they are not confidence intervals.\n'
                 'Nominal associations do not establish multiplicity-adjusted evidence or causation.',
                 fontsize=10, linespacing=1.8)
        return fig


def save_figure(fig, output_dir, stem, dpi=300):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    with plt.rc_context(STYLE):
        for extension in ['png', 'pdf', 'svg']:
            path = output_dir / f'{stem}.{extension}'
            fig.savefig(path, dpi=dpi, facecolor='white', bbox_inches='tight', pad_inches=.18)
            paths.append(path)
    return paths


def export_supporting_files(workbook, studies, alkan, garcia, output_dir):
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    captions = output_dir / 'Figure_captions.txt'
    captions.write_text(CAPTION_2 + '\n\n' + CAPTION_3 + '\n', encoding='utf-8')
    samples_csv = output_dir / 'Figure_2_plotted_data.csv'
    # Export only fields relevant to the bars; preserve source locations.
    columns = ['study_id', 'study_label', 'total_reported_n', 'microbiome_n',
               'sdf_n', 'sdf_n_status', 'doi', 'source', 'excel_row']
    studies[columns].to_csv(samples_csv, index=False, encoding='utf-8-sig', na_rep='')
    correlations_csv = output_dir / 'Figure_3_plotted_data.csv'
    pd.concat([alkan.assign(panel='A'), garcia.assign(panel='B')], ignore_index=True).to_csv(
        correlations_csv, index=False, encoding='utf-8-sig', na_rep='')
    metadata = output_dir / 'Run_information.json'
    metadata.write_text(json.dumps({
        'input_workbook': Path(workbook).name,
        'input_sha256': hashlib.sha256(Path(workbook).read_bytes()).hexdigest(),
        'python': platform.python_version(), 'pandas': pd.__version__,
        'numpy': np.__version__, 'matplotlib': matplotlib.__version__,
        'openpyxl': openpyxl.__version__,
        'plotted_rows': {'Samples': len(studies), 'Alkan': len(alkan), 'Garcia_Segura': len(garcia)},
        'analysis': 'Descriptive figures only; no pooled effects or confidence intervals.',
    }, indent=2), encoding='utf-8')
    return [captions, samples_csv, correlations_csv, metadata]


def make_zip(paths, zip_path):
    """Include only this run's explicit outputs, avoiding stale files."""
    zip_path = Path(zip_path)
    with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as archive:
        for path in paths:
            archive.write(path, arcname=Path(path).name)
    return zip_path


def upload_workbook():
    from google.colab import files
    print('Choose your Study_Data_Extraction Excel workbook (.xlsx).')
    uploaded = files.upload()
    choices = [name for name in uploaded if name.lower().endswith('.xlsx')]
    if len(choices) != 1:
        raise ValueError('Please rerun this cell and upload exactly one .xlsx workbook.')
    return Path(choices[0])


def main(workbook=None, output_dir='semen_microbiome_figures', dpi=300, show=True):
    in_colab = 'google.colab' in sys.modules
    if workbook is None:
        if not in_colab:
            raise ValueError('Supply the path to your extraction workbook.')
        workbook = upload_workbook()
    studies, alkan, garcia = read_workbook(workbook)
    paths = []
    for fig, stem in [
        (draw_coverage(studies), 'Figure_2_Sample_coverage'),
        (draw_correlations(alkan, garcia), 'Figure_3_Reported_correlations'),
    ]:
        paths.extend(save_figure(fig, output_dir, stem, dpi))
        if show:
            plt.show()
        plt.close(fig)
    paths.extend(export_supporting_files(workbook, studies, alkan, garcia, output_dir))
    archive = make_zip(paths, Path(output_dir).with_suffix('.zip'))
    print(f'Created Figures 2 and 3, captions and plotted data: {archive}')
    if in_colab:
        from google.colab import files
        files.download(str(archive))
    return archive


if __name__ == '__main__':
    if 'google.colab' in sys.modules:
        main()
    else:
        import argparse
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('workbook', help='Path to the extraction .xlsx workbook')
        parser.add_argument('--output-dir', default='semen_microbiome_figures')
        parser.add_argument('--dpi', type=int, default=300)
        args = parser.parse_args()
        main(args.workbook, args.output_dir, args.dpi, show=False)
