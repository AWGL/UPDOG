#!/usr/bin/env python

import io
import pandas as pd
from pyvariantfilter.variant import Variant
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import scipy

sns.set()

COUNTER_KEYS = [
    'variant_count',
    'matches_paternal_uniparental_ambiguous_count',
    'matches_maternal_uniparental_ambiguous_count',
    'matches_paternal_uniparental_isodisomy_count',
    'matches_maternal_uniparental_isodisomy_count',
    'alleles_identical_to_dad_count',
    'alleles_identical_to_mum_count',
    'is_homozygous_count',
    'is_biparental_count',
]

def extract_genotype_data(sample_genotype_data, ref, alt, chrom):
	"""
	Extract and clean GT/AD/GQ/DP/AF for one sample at one variant
	"""
	
	ref_and_alt = [ref, alt]
	gts = ['.' if allele is None else ref_and_alt[allele]
    	for allele in sample_genotype_data['GT']]

	if len(gts) == 1 and chrom.lstrip('chr') in ('X', 'Y'):
		gts.append(gts[0])

	if gts == ['.', '.'] or (len(sample_genotype_data['AD']) == 1 and sample_genotype_data['AD'][0] is None):
		ads = [0, 0]
	else:
		assert len(sample_genotype_data['AD']) == 2
		ads = [0 if ad is None else ad for ad in sample_genotype_data['AD']]

	gq = sample_genotype_data['GQ'] or 0
	dp = sample_genotype_data['DP'] or 0
	af_raw = sample_genotype_data['AF']
	af = 0 if af_raw is None else af_raw[0]

	return gts, ads, gq, dp, af

def variant_passes_prefilters(alt_alleles, alt, quality, min_qual):
    """
    True if the variant should be processed further
    """
    if len(alt_alleles) != 1:
        return False
    if alt == '*':
        return False
    if quality <= min_qual:
        return False
    return True

def has_crossed_block_boundary(pos, last_block, block_size):
    return pos > last_block + block_size

def init_counters():
    return {key: 0 for key in COUNTER_KEYS}

def accumulate_counters(counters, flags):
    """
    flags: dict of the same keys (minus 'variant_count') mapped to bool
    """
    counters['variant_count'] += 1
    for key, matched in flags.items():
        if matched:
            counters[key] += 1
    return counters

def finalize_block(counters, block_count, block_size):
    return {'end': block_count * block_size, **counters}

def get_variant_match_flags(variant, proband_id, min_gq, min_dp):
    return {
        'matches_paternal_uniparental_ambiguous_count': variant.matches_paternal_uniparental_ambiguous(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'matches_maternal_uniparental_ambiguous_count': variant.matches_maternal_uniparental_ambiguous(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'matches_paternal_uniparental_isodisomy_count': variant.matches_paternal_uniparental_isodisomy(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'matches_maternal_uniparental_isodisomy_count': variant.matches_maternal_uniparental_isodisomy(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'alleles_identical_to_dad_count': variant.alleles_identical_to_dad(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'alleles_identical_to_mum_count': variant.alleles_identical_to_mum(min_parental_gq=min_gq, min_parental_depth=min_dp),
        'is_homozygous_count': variant.is_hom_alt(proband_id),
        'is_biparental_count': variant.is_biparental_inheritance(min_parental_gq=min_gq, min_parental_depth=min_dp),
    }

def calculate_upd_metrics_per_chromosome(bcf_in, chromosome_to_analyze, family, block_size, min_dp, min_gq, min_qual, proband_id):
	"""
	Main function to calculate upd metrics
	"""
	family_member_ids = family.get_all_family_member_ids()
	counters = init_counters()
	block_dict = {}
	variant_list = []
	last_block = 0
	block_count = 1

	for rec in bcf_in.fetch(chromosome_to_analyze):
		if not variant_passes_prefilters(rec.alts, rec.alts[0] if rec.alts else None, rec.qual, min_qual):
			continue

		alt = rec.alts[0]
		new_variant = Variant(chrom=rec.chrom.strip('chr'), pos=rec.pos, ref=rec.ref, alt=alt,
                               filter_status=rec.filter.keys(), quality=rec.qual)
		new_variant.add_family(family)

		sample_data = {}
		for family_member_id in family_member_ids:
			gts, ads, gq, dp, af = extract_genotype_data(rec.samples[family_member_id], rec.ref, alt, rec.chrom)
			new_variant.add_genotype(family_member_id, gts, ads, gq, dp)
			sample_data[family_member_id] = (ads, gq, dp, af)

		if not (new_variant.is_snp() and new_variant.all_samples_pass_genotype_quality(min_dp, min_gq)):
			continue
		
		#Get metrics just for proband
		p_ads, p_gq, p_dp, p_af = sample_data[proband_id]
		variant_list.append([rec.chrom, rec.pos, p_ads[0], p_ads[1], p_af, p_dp, p_gq, rec.qual])
		flags = get_variant_match_flags(new_variant, proband_id, min_gq, min_dp)
		accumulate_counters(counters, flags)

		if has_crossed_block_boundary(new_variant.pos, last_block, block_size):
			block_dict[block_count] = finalize_block(counters, block_count, block_size)
			counters = init_counters()
			last_block += block_size
			block_count += 1

	return block_dict, variant_list

def replace_with_na(df, column, min_count):

	if int(df['variant_count']) < int(min_count):
		
		return np.nan
		
	else:
		
		return df[column]

def create_ax_for_plotting(chromosome, df, block_size, output):

	chrom_prop_df = df[df['chrom'] ==chromosome]

	plot_min = min(chrom_prop_df['end'])
	plot_max = max(chrom_prop_df['end'])

	# we want to scale the X axis labels so that it looks pretty
	chrom_size = plot_max

	if chrom_size > 100000000:
		
		scale_factor = 10
		
	else:
		
		scale_factor = 3

	xticks_chrom = np.arange(plot_min, plot_max, block_size*scale_factor)
	x_ticks_labels = [int(x/1000000) for x in xticks_chrom]

	# plot and format axis
	ax = chrom_prop_df.plot(x='end', figsize=(20,5))
	ax.set_xticks(xticks_chrom)
	ax.set_xlim([plot_min-1000000, plot_max+1000000])
	ax.set_ylim([-0.02, 1.05])
	ax.set_xticklabels(x_ticks_labels)
	ax.set_xlabel(f'Chromosome {chromosome} Position (Mb)')
	ax.set_ylabel('Proportion of Variants')
	ax.legend(loc='upper right')

	plt.savefig(output, format='png')

	plt.close('all')

def plot_variants(chromosome, df, output, block_size, sample_fraction):
	"""
	Plots allele frequency vs. position for a given chromosome with optional downsampling.

	Parameters:
	chromosome (str): The chromosome to plot.
	df (pd.DataFrame): DataFrame containing the variant data with 'chrom', 'pos', and 'af' columns.
	output (str): Path to save the plot as a PNG file.
	block_size (int): Block size for xticks calculation.
	sample_fraction (float): Fraction of data to sample for plotting (default: 0.1).
	"""
	# Filter data for the specific chromosome

	chrom_prop_df = df[df['chrom'] == chromosome]

	# Downsample the data
	if sample_fraction < 1.0:
		chrom_prop_df = chrom_prop_df.sample(frac=sample_fraction, random_state=42)

	# Calculate plot range and scaling
	plot_min = chrom_prop_df['pos'].min()
	plot_max = chrom_prop_df['pos'].max()

	chrom_size = plot_max

	scale_factor = 10 if chrom_size > 100000000 else 3

	xticks_chrom = np.arange(plot_min, plot_max, block_size * scale_factor)
	x_ticks_labels = [int(x / 1000000) for x in xticks_chrom]

	# Plot BAF plot with a smaller depth plot sharing the same X-axis directly underneath
	fig, (ax, ax_depth) = plt.subplots(2, 1, figsize=(20, 6), sharex=True, gridspec_kw={'height_ratios': [3, 1]})
	
	#Set limits for depth plot
	mean_dp = chrom_prop_df['dp'].mean()
	if pd.notna(mean_dp) and mean_dp > 0:
		depth_y_lim = mean_dp * 3
	else:
		depth_y_lim = 150

	#Do the plotting
	ax.scatter(x='pos', y='af', data=chrom_prop_df, s=2, alpha=0.5)
	ax.set_xlim([plot_min - 1000000, plot_max + 1000000])
	ax.set_ylim([-0.02, 1.05])
	ax.set_ylabel('Beta Allele Frequency')
	ax.tick_params(labelbottom=False)

	ax_depth.scatter(x='pos', y='dp', data=chrom_prop_df, s=2, alpha=0.5, color='grey')
	ax_depth.set_xticks(xticks_chrom)
	ax_depth.set_xticklabels(x_ticks_labels)
	ax_depth.set_xlabel(f'Chromosome {chromosome} Position (Mb)')
	ax_depth.set_ylabel('Depth')
	ax_depth.set_ylim(bottom=0, top=depth_y_lim)

	fig.subplots_adjust(hspace=0.05)

	# Save and close plot
	plt.savefig(output, format='png')
	plt.close('all')


def is_significant(df, expected, key):
	
	if df[key] == None:
		
		return None
	
	mendel_errors = df[key]
	
	variant_count = df['variant_count']
	
	
	return scipy.stats.binomtest(mendel_errors, variant_count, expected, alternative='greater').pvalue


def merge_contiguous_blocks(df, block_size, analysis, analysis2):

	if df.shape[0] == 0:

		return []

	contiguous_blocks = []

	last_chr = 'NA'
	last_start = None
	last_end = None

	current_block_start = None
	current_p_values = []
	current_proportions = []

	for row in df.itertuples():
		
		current_chr = row.chrom
		current_end = row.end
		current_start = current_end - block_size

		if current_chr != last_chr:
			
			if last_chr != 'NA':
				
				contiguous_blocks.append([last_chr, current_block_start, last_end, np.mean(current_p_values), np.mean(current_proportions)])
				
				current_p_values = []
				
			current_block_start = current_start
			last_chr = current_chr
			
		elif current_start != (last_start+block_size):
			
			contiguous_blocks.append([current_chr, current_block_start, last_end, np.mean(current_p_values), np.mean(current_proportions)])
			current_block_start = current_start
			current_p_values  =[]
			
		
		last_chr = current_chr
		last_start = current_start
		last_end = current_end
		current_p_values.append(row.__getattribute__(analysis))
		current_proportions.append(row.__getattribute__(analysis2))
			
	contiguous_blocks.append([last_chr, current_block_start, last_end, np.mean(current_p_values), np.mean(current_proportions)])

	return contiguous_blocks

def apply_filters(df, min_blocks, min_proportion, block_size):

	filters = []

	proportion = df['mean_proportion_me']
	start = df['start']
	end = df['end']

	if proportion < min_proportion:

		filters.append('low_proportion_me')

	block_length = end - start

	n_blocks = block_length / block_size

	if n_blocks < min_blocks:

		filters.append('min_block_length')

	if len(filters) == 0:

		return 'pass'

	else:

		return ';'.join(filters)

def get_genome_build(header):
	"""
	Get the genome build from the header
	"""
	for line in io.StringIO(header).readlines():

		if '##reference' in line:

			if 'hs37d5' in line or 'human_g1k_v37' in line or 'GRCh37' in line or 'hg19' in line:

				return 37

			elif 'Homo_sapiens_assembly38' in line or 'GRCh38' in line or 'hg38' in line:

				return 38

	raise Exception('Invalid reference genome in VCF header')

def get_info_from_ped(ped, proband_id):
	"""
	Parse out all information needed from ped file
	
	:param ped: ped file
	:param proband_id: sample identifier for proband
	"""
	# read ped into df
	ped_df = pd.read_csv(ped, sep='\t', names=['family_id', 'sample_id', 'paternal_id', 'maternal_id', 'sex', 'affected'], dtype={'sex': str,
																																'affected': str,
																																'paternal_id': str,
																																'maternal_id': str,
																																'family_id': str,
																																})
	# filter by proband
	filtered_ped = ped_df[ped_df['sample_id']==proband_id]

	# get mum and dad ids
	dad = filtered_ped['paternal_id'].iloc[0]
	mum = filtered_ped['maternal_id'].iloc[0]
	sex = filtered_ped['sex'].iloc[0]
	affected = filtered_ped['affected'].iloc[0]

	# get family id
	family_id = filtered_ped['family_id'].iloc[0]

	return dad, mum, sex, affected, family_id

def singleton_status(dad, mum):
	"""
	Determine if singleton if either mum or dad, or both, missing
	
	:param dad: Dad sample ID
	:param mum: Mum sample ID
	"""
	# Record as singleton if we don't have mum and dad
	if dad == '0' or mum == '0':

		return True

	else:

		return False
	
def set_chromosomes(just_one_chromosome, chromosome, sex, genome_build):
	"""
	Set chromosomes to analyse
	"""
	# if set to a single chromosome
	if just_one_chromosome == True:

		if sex == '1' and chromosome == 'X':

			raise ValueError('Chromosome cannot be X if proband is male.')

		else:

			chromosomes_to_analyze = [chromosome]

	else:

		chromosomes_to_analyze = ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X']

		if sex == '2':

			print ('Proband is Female - will analyse chromosomes 1-22 and X')

		elif sex == '1':

			print ('Proband is Male - will analyse chromosomes 1-22 only.')
			chromosomes_to_analyze.pop()

	if genome_build == 38:
			
		print('Genome build is 38')
		return ['chr' + c for c in chromosomes_to_analyze]
	
	elif genome_build == 37:
		print('Genome build is 37')
		return chromosomes_to_analyze
	
def check_affected_and_sex(affected, sex):
	"""
	Simple function to make sure that affected and sex are set - triggers sys.exit()
	"""

	if affected != '2':

		return True
	
	if sex == '0':

		return True
	
	else:

		return False

