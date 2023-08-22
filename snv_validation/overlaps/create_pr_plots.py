import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import argparse
import os

# Genotype categories
GT_CATS = ['./.', '0/0', '0/1', '1/1']

# Columns in input file
COLS = [
    'chr',
    'pos',
    'rs',
    'dbsnp_common',
    'dbsnp_max_maf',
    'gt_bisc',
    'gt_giab',
    'gq_bisc',
    'gq_giab',
    'af_bisc',
    'af_giab',
    'dp_bisc',
    'dp_giab',
]

# Minimum genotype quality for dbSNP overlap
DBSNP_MIN_SCORE = 15

# Number of "bins" or genotype qualities
N_BINS=257

# Points to mark
POINTS = [15, 60]

def parse_cli():
    """Parse command line arguments.

    Inputs -
        None
    Returns -
        argparser.ArgumentParser.parse_args() object
    """
    parser = argparse.ArgumentParser(description='Create PR curves for BISCUIT SNPs')

    parser.add_argument('input', help='Output file from preprocess_files.sh')
    parser.add_argument(
        '-d', '--dbsnp', type=int, choices=[1,2,3,4], default=4,
        help='1: ignore dbSNP info, 2: remove dbSNP SNPs, 3: only use dbSNP SNPs, 4: use dbSNP info for filtering [default]'
    )
    parser.add_argument(
        '-g', '--giab', action='store_true',
        help='Keep only SNPs with GIAB GQ == 99 or GIAB GQ == NA'
    )

    return parser.parse_args()

def pr_curve_points(df, masks):
    """Create x and y points for the PR curves.

    Inputs -
        df    - DataFrame to extract data from
        masks - masks on the DataFrame with TP/FP/FN results
    Returns -
        tuple, (x points, y points)
    """
    tp_counts = [0 for _ in range(N_BINS)]
    fp_counts = [0 for _ in range(N_BINS)]
    fn_counts = [0 for _ in range(N_BINS)]

    # Iterate backwards to make it cumulative
    x = []
    y = []
    for bin in range(N_BINS-1, -1, -1):
        tp_counts[bin] = np.nansum(df[masks[0]]['gq_bisc'] == bin)
        fp_counts[bin] = np.nansum(df[(masks[1]) | (masks[2])]['gq_bisc'] == bin)
        fn_counts[bin] = np.nansum(df[masks[3]]['gq_bisc'] == bin) # shouldn't include any uncallable as they have score NA

        # On the first bin, include FN_uncallable
        if bin == N_BINS-1:
            fn_counts[bin] += sum(masks[4])
        else:
            tp_counts[bin] += tp_counts[bin+1]
            fp_counts[bin] += fp_counts[bin+1]
            fn_counts[bin] += fn_counts[bin+1]

        if tp_counts[bin] + fn_counts[bin] == 0:
            x.append(np.nan)
        else:
            x.append(tp_counts[bin] / (tp_counts[bin] + fn_counts[bin]))
        if tp_counts[bin] + fp_counts[bin] == 0:
            y.append(np.nan)
        else:
            y.append(tp_counts[bin] / (tp_counts[bin] + fp_counts[bin]))

    return (x, y)

def create_plot(orig_x, orig_y, orig_l, filt_x, filt_y, filt_l, title, outname):
    """Create PR plot.

    Inputs -
        orig_x  - x-points for original data
        orig_y  - y-points for original data
        orig_l  - label for original data
        filt_x  - x-points for filtered data
        filt_y  - y-points for filtered data
        filt_l  - label for filtered data
        title   - figure title
        outname - output file name
    Returns -
        None
    """
    fig = plt.figure(figsize=(7, 7))

    plt.plot(orig_x, orig_y, 'k-', label=orig_l)
    plt.plot(filt_x, filt_y, 'r-', label=filt_l)

    for idx in range(N_BINS):
        tmp = N_BINS - idx - 1
        xo_shift = -0.03 if tmp < 40 else 0
        yo_shift = 0 if tmp < 40 else  0.03
        xf_shift =  0.03 if tmp < 40 else 0
        yf_shift = 0 
        if tmp in POINTS and orig_x[idx] != 0:
            plt.plot(orig_x[idx], orig_y[idx], 'ko')
            plt.plot(filt_x[idx], filt_y[idx], 'ro')

        if tmp == 15:
            print(f'15: {orig_y[idx]} {filt_y[idx]}')
            plt.annotate(tmp, (orig_x[idx]-0.07, orig_y[idx]-0.01), fontsize=16, color='black')
            plt.annotate(tmp, (filt_x[idx]+0.01, filt_y[idx]-0.01), fontsize=16, color='red')
        elif tmp == 60:
            print(f'60: {orig_y[idx]} {filt_y[idx]}')
            plt.annotate(tmp, (orig_x[idx]-0.01, orig_y[idx]+0.005), fontsize=16, color='black')
            plt.annotate(tmp, (filt_x[idx]+0.01, filt_y[idx]+0.005), fontsize=16, color='red')

    plt.legend(loc='lower left', fontsize=20)
    plt.grid()

    plt.title(title, fontsize=24)
    plt.xlabel('Recall'   , fontsize=20)
    plt.ylabel('Precision', fontsize=20)

    plt.xlim(-0.03, 1.03)
    plt.ylim( 0.69, 1.01)

    plt.xticks(
        [i for i in np.arange(0.0, 1.01, 0.1)],
        [str(round(i,1)) for i in np.arange(0.0, 1.01, 0.1)],
        fontsize=18, rotation=90
    )
    plt.yticks(
        [i for i in np.arange(0.7, 1.01, 0.1)],
        [str(round(i,1)) for i in np.arange(0.7, 1.01, 0.1)],
        fontsize=18
    )

    plt.savefig(outname, bbox_inches='tight')
    plt.close('all')

    return None

def main():
    """Process files.

    Inputs -
        None
    Returns -
        None
    """
    args = parse_cli()

    name = os.path.basename(args.input).replace('.processed.tsv', '')

    # Read data
    df_orig = pd.read_csv(args.input, sep='\t', header=None, names=COLS, na_values=['.', 'NA', 'N/A'])

    # Remove any genotypes that we don't want (1/2, 2/3, etc.) - these are a small fraction of the results
    df_orig = df_orig[(df_orig['gt_bisc'].isin(GT_CATS)) & (df_orig['gt_giab'].isin(GT_CATS))]

    # If requested by user, apply filtering for GIAB GQ values
    if args.giab:
        df_orig = df_orig[(df_orig['gq_giab'] == 99) | (df_orig['gq_giab'].isna())]

    # Apply filters based on dbSNP results
    df_filt = df_orig.copy(deep=True)
    if args.dbsnp == 4:
        dbsnp    = (~df_filt['rs'].isna()) & (~df_filt['dbsnp_common'].isna()) & (df_filt['dbsnp_max_maf'] > 0.05)
        name     = f'{name}.dbsnp{DBSNP_MIN_SCORE}'
        meetsmin = (df_filt['gq_bisc'] >= DBSNP_MIN_SCORE)

        df_filt.loc[(dbsnp & meetsmin), 'gq_bisc'] = 255
    elif args.dbsnp == 3:
        dbsnp = (~df_filt['rs'].isna()) & (~df_filt['dbsnp_common'].isna()) & (df_filt['dbsnp_max_maf'] > 0.05)
        name  = f'{name}.dbsnponly'
        df_filt    = df_filt[dbsnp]
    elif args.dbsnp == 2:
        dbsnp = (~df_filt['rs'].isna())
        name  = f'{name}.notdbsnp'
        df_filt    = df_filt[~dbsnp]

    dfs = [df_orig, df_filt] # unfiltered results, dbSNP prior results

    # Set up masks for true positive / true negative categories
    masks = [[] for _ in range(4)] # 2 sets (het, hom) for each data frame (original, filtered)
    for i in range(len(dfs)):
        df = dfs[i]

        # Should 1/1 be included in FN for GIAB?
        TP_het     = (df['gt_bisc'] == '0/1') & (df['gt_giab'] == '0/1')
        FP_het_ref = (df['gt_bisc'] == '0/1') & ( (df['gt_giab'] == '0/0') | (df['gt_giab'] == './.') )
        FP_hom     = (df['gt_bisc'] == '0/1') & (df['gt_giab'] == '1/1')
        FN_het     = ( (df['gt_bisc'] == '0/0') | (df['gt_bisc'] == '1/1') ) & (df['gt_giab'] == '0/1')
        FN_het_uncallable = (df['gt_bisc'] == './.') & (df['gt_giab'] == '0/1')

        masks[2*i] = [TP_het, FP_het_ref, FP_hom, FN_het, FN_het_uncallable]
        
        # Should 0/1 be included in FN for GIAB?
        TP_hom     = (df['gt_bisc'] == '1/1') & (df['gt_giab'] == '1/1')
        FP_hom_ref = (df['gt_bisc'] == '1/1') & ( (df['gt_giab'] == '0/0') | (df['gt_giab'] == './.') )
        FP_het     = (df['gt_bisc'] == '1/1') & (df['gt_giab'] == '0/1')
        FN_hom     = ( (df['gt_bisc'] == '0/0') | (df['gt_bisc'] == '0/1') ) & (df['gt_giab'] == '1/1')
        FN_hom_uncallable = (df['gt_bisc'] == './.') & (df['gt_giab'] == '1/1')

        masks[2*i+1] = [TP_hom, FP_hom_ref, FP_het, FN_hom, FN_hom_uncallable]

    # Precision-Recall plots
    if True:
        orig_x_het, orig_y_het = pr_curve_points(dfs[0], masks[0])
        orig_x_hom, orig_y_hom = pr_curve_points(dfs[0], masks[1])
        filt_x_het, filt_y_het = pr_curve_points(dfs[1], masks[2])
        filt_x_hom, filt_y_hom = pr_curve_points(dfs[1], masks[3])

        create_plot(
            orig_x_het, orig_y_het, 'GQ >= n',
            filt_x_het, filt_y_het, 'GQ >= 15 (dbSNP+)\nGQ >= n (dbSNP-)',
            'chr11p15 Heterozygous SNPs',
            'precision_recall_het.pdf'
        )
        create_plot(
            orig_x_hom, orig_y_hom, 'GQ >= n',
            filt_x_hom, filt_y_hom, 'GQ >= 15 (dbSNP+)\nGQ >= n (dbSNP-)',
            'chr11p15 Homozygous SNPs',
            'precision_recall_hom.pdf'
        )

    return None

if __name__ == '__main__':
    main()
