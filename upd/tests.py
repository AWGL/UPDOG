import unittest
from unittest.mock import MagicMock

from utility_funcs import *

class test_utility_funcs(unittest.TestCase):
    """
    Tests for utility functions
    """

    def test_genome_build(self):
        """
        Test assignment of genome build
        """
        # 'hs37d5'
        header = "###reference=hs37d5"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 37)

        # 'human_g1k_v37' 
        header = "###reference=human_g1k_v37"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 37)

        # 'GRCh37' 
        header = "###reference=GRCh37"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 37)

        # 'hg19'
        header = "###reference=hg19"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 37)

        #'Homo_sapiens_assembly38' 
        header = "###reference=Homo_sapiens_assembly38"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 38)

        # 'GRCh38'
        header = "###reference=GRCh38"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 38)

        # 'hg38'
        header = "###reference=hg38"
        genome_build = get_genome_build(header)
        self.assertEqual(genome_build, 38)

        #Incorrect header
        header = "###reference=nonsense"
        with self.assertRaises(Exception) as ctx:
            get_genome_build(header)
        self.assertEqual(str(ctx.exception), "Invalid reference genome in VCF header")
      
        # No header
        header = ""
        with self.assertRaises(Exception) as ctx:
            get_genome_build(header)
        self.assertEqual(str(ctx.exception), "Invalid reference genome in VCF header")

        # Header, just no reference
        header = "some other string"
        with self.assertRaises(Exception) as ctx:
            get_genome_build(header)
        self.assertEqual(str(ctx.exception), "Invalid reference genome in VCF header")

    def test_parse_ped(self):
        """
        Test for parsing ped to get information
        """

        #Test simple ped with single family
        ped = io.StringIO(
            "FAM001\tPROBAND1\tDAD1\tMUM1\t1\t2\n"
            "FAM001\tDAD1\t0\t0\t1\t1\n"
            "FAM001\tMUM1\t0\t0\t2\t1\n"
        )

        dad, mum, sex, affected, family_id = get_info_from_ped(ped, "PROBAND1")

        self.assertEqual(dad, "DAD1")
        self.assertEqual(mum, "MUM1")
        self.assertEqual(sex, "1")
        self.assertEqual(affected, "2")
        self.assertEqual(family_id, "FAM001")

        #More than one family in PED, get the right one
        ped = io.StringIO(
            "FAM001\tPROBAND1\tDAD1\tMUM1\t1\t2\n"
            "FAM001\tDAD1\t0\t0\t1\t1\n"
            "FAM001\tMUM1\t0\t0\t2\t1\n"
            "FAM002\tPROBAND2\tDAD2\tMUM2\t2\t2\n"
            "FAM002\tDAD2\t0\t0\t1\t1\n"
            "FAM002\tMUM2\t0\t0\t2\t1\n"
            "FAM003\tPROBAND3\t0\t0\t1\t2\n"
        )

        dad, mum, sex, affected, family_id = get_info_from_ped(ped, "PROBAND1")

        self.assertEqual(dad, "DAD1")
        self.assertEqual(mum, "MUM1")
        self.assertEqual(sex, "1")
        self.assertEqual(affected, "2")
        self.assertEqual(family_id, "FAM001")

        #Check that singleton info is pulled correctly too
        ped = io.StringIO(
            "FAM001\tPROBAND1\tDAD1\tMUM1\t1\t2\n"
            "FAM001\tDAD1\t0\t0\t1\t1\n"
            "FAM001\tMUM1\t0\t0\t2\t1\n"
            "FAM002\tPROBAND2\tDAD2\tMUM2\t2\t2\n"
            "FAM002\tDAD2\t0\t0\t1\t1\n"
            "FAM002\tMUM2\t0\t0\t2\t1\n"
            "FAM003\tPROBAND3\t0\t0\t1\t2\n"
        )

        dad, mum, sex, affected, family_id = get_info_from_ped(ped, "PROBAND3")

        self.assertEqual(dad, "0")
        self.assertEqual(mum, "0")
        self.assertEqual(sex, "1")
        self.assertEqual(affected, "2")
        self.assertEqual(family_id, "FAM003")

    def test_singleton_status(self):
        """
        Test function that determines if singleton based on parental IDs
        """

        #Mum and dad - false
        mum = "MUM1"
        dad = "DAD1"

        singleton = singleton_status(dad, mum)
        self.assertFalse(singleton)

        #Just mum - true
        dad = "0"
        
        singleton = singleton_status(dad, mum)
        self.assertTrue(singleton)

        #Just dad - true
        dad = "DAD1"
        mum = "0"

        singleton = singleton_status(dad, mum)
        self.assertTrue(singleton)

        #Neither - true
        dad = "0"

        singleton = singleton_status(dad, mum)
        self.assertTrue(singleton)

    def test_set_chromosomes(self):
        """
        Test assignment of chromosomes to analyse
        """

        # All chromosomes, male, b38
        just_one_chromosome = False
        chromosome = None
        sex = "1"
        genome_build = 38
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22'])

        # All chromosomes, male, b37
        genome_build = 37
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22'])

        # All chromosomes, female, b38
        sex = "2"
        genome_build = 38
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX'])

        # All chromosomes, female, b37
        genome_build = 37
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '20', '21', '22', 'X'])

        # Single chromosome, non-sex, male, b38
        just_one_chromosome = True
        chromosome = "3"
        sex = "1"
        genome_build = 38
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['chr3'])

        # Single chromosome, non-sex, male, b37
        genome_build = 37
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['3'])

        # Single chromosome, non-sex, female, b38
        sex = "2"
        genome_build = 38
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['chr3'])

        # Single chromosome, non-sex, female, b37
        genome_build = 37
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['3'])
        
        # Single chromosome, X, male, b38
        chromosome = "X"
        sex = "1"
        genome_build = 38
        with self.assertRaises(Exception) as ctx:
            set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(str(ctx.exception), "Chromosome cannot be X if proband is male.")

        # Single chromosome, X, male, b37
        genome_build = 37
        with self.assertRaises(Exception) as ctx:
            set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(str(ctx.exception), "Chromosome cannot be X if proband is male.")

        # Single chromosome, X, female, b38
        sex = "2"
        genome_build = 38
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['chrX'])

        # Single chromosome, X, female, b37
        genome_build = 37
        chromosomes_to_analyze = set_chromosomes(just_one_chromosome, chromosome, sex, genome_build)
        self.assertEqual(chromosomes_to_analyze, ['X'])

    def test_extract_genotype_data(self):
        """
        Testing the extraction of genotype data
        """
        
        #No values
        sample = {'GT': (None, None), 'AD': (None,), 'GQ': None, 'DP': None, 'AF': None}
        gts, ads, gq, dp, af = extract_genotype_data(sample, ref='A', alt='T', chrom='1')
        self.assertEqual(gts, ['.', '.'])
        self.assertEqual(ads, [0, 0])
        self.assertEqual(gq, 0)
        self.assertEqual(dp, 0)
        self.assertEqual(af, 0)

        #gts is none and something else
        sample = {'GT': (1, None), 'AD': (None,), 'GQ': None, 'DP': None, 'AF': None}
        gts, ads, gq, dp, af = extract_genotype_data(sample, ref='A', alt='T', chrom='1')
        self.assertEqual(gts, ['T', '.'])

        #Hemizygous (on X or Y), the genotype is doubled
        sample = {'GT': (1,), 'AD': (3, 7), 'GQ': 40, 'DP': 10, 'AF': [0.7]}
        gts, ads, gq, dp, af = extract_genotype_data(sample, ref='A', alt='T', chrom='chrX')
        self.assertEqual(gts, ['T', 'T'])  # doubled for hemizygous

    def test_pass_prefilter(self):
        """
        Testing the filtering of variants not included
        """

        #All good
        alt_alleles = ("T")
        alt = ("T")
        quality = 40
        min_qual= 5
        pass_var = variant_passes_prefilters(alt_alleles, alt, quality, min_qual)
        self.assertTrue(pass_var)

        #Too many alts
        alt_alleles = ("T","A")
        pass_var = variant_passes_prefilters(alt_alleles, alt, quality, min_qual)
        self.assertFalse(pass_var)

        #Alt is *
        alt_alleles = ("*")
        alt = ("*")
        pass_var = variant_passes_prefilters(alt_alleles, alt, quality, min_qual)
        self.assertFalse(pass_var)

        #Quality is below cut off
        alt_alleles = ("T")
        alt = ("T")
        quality = 4
        min_qual= 5
        pass_var = variant_passes_prefilters(alt_alleles, alt, quality, min_qual)
        self.assertFalse(pass_var)

    def test_cross_block_boundary(self):
        """
        Testing the check for crossing a block boundary
        """

        pos = 100
        last_block = 5
        block_size = 5
        check = has_crossed_block_boundary(pos, last_block, block_size)
        self.assertTrue(check)

        pos = 1
        check = has_crossed_block_boundary(pos, last_block, block_size)
        self.assertFalse(check)

    def test_accumulate_counters(self):
        """
        Testing the accumulating counters
        """
        counters = init_counters()
        flags = {
            'matches_paternal_uniparental_ambiguous_count': True,
            'matches_maternal_uniparental_ambiguous_count': False,
            'matches_paternal_uniparental_isodisomy_count': False,
            'matches_maternal_uniparental_isodisomy_count': False,
            'alleles_identical_to_dad_count': True,
            'alleles_identical_to_mum_count': False,
            'is_homozygous_count': False,
            'is_biparental_count': False,
        }
        accumulate_counters(counters, flags)
        self.assertEqual(counters['variant_count'], 1)
        self.assertEqual(counters['matches_paternal_uniparental_ambiguous_count'], 1)
        self.assertEqual(counters['alleles_identical_to_dad_count'], 1)
        self.assertEqual(counters['matches_maternal_uniparental_ambiguous_count'], 0)

        # Second call on the same dict - confirm it accumulates, not resets
        flags2 = {
            'matches_paternal_uniparental_ambiguous_count': False,
            'matches_maternal_uniparental_ambiguous_count': True,
            'matches_paternal_uniparental_isodisomy_count': False,
            'matches_maternal_uniparental_isodisomy_count': False,
            'alleles_identical_to_dad_count': True,
            'alleles_identical_to_mum_count': False,
            'is_homozygous_count': False,
            'is_biparental_count': False,
        }
        accumulate_counters(counters, flags2)
        self.assertEqual(counters['variant_count'], 2)
        self.assertEqual(counters['matches_paternal_uniparental_ambiguous_count'], 1)  # unchanged from call 1
        self.assertEqual(counters['matches_maternal_uniparental_ambiguous_count'], 1)  # newly incremented
        self.assertEqual(counters['alleles_identical_to_dad_count'], 2)  # incremented both times

    def test_finalize_block(self):
    
        counters = {'variant_count': 5, 'is_homozygous_count': 2}
        result = finalize_block(counters, block_count=3, block_size=1000000)
        self.assertEqual(result, {'end': 3000000, 'variant_count': 5, 'is_homozygous_count': 2})

    def test_replace_with_na(self):

        df = {'variant_count': 5, 'mean_proportion_me': 0.3}
        self.assertTrue(np.isnan(replace_with_na(df, 'mean_proportion_me', min_count=10)))
        self.assertEqual(replace_with_na(df, 'mean_proportion_me', min_count=5), 0.3) 

    def test_is_significant(self):

        df = {'mendel_errors': None, 'variant_count': 100}
        self.assertIsNone(is_significant(df, expected=0.01, key='mendel_errors'))

        df = {'mendel_errors': 50, 'variant_count': 100}
        pvalue = is_significant(df, expected=0.01, key='mendel_errors')
        self.assertLess(pvalue, 0.001)  

    def test_apply_filters(self):

        # passes both checks
        df = {'mean_proportion_me': 0.5, 'start': 0, 'end': 10000000}
        self.assertEqual(apply_filters(df, min_blocks=5, min_proportion=0.01, block_size=1000000), 'pass')

        # low proportion only
        df = {'mean_proportion_me': 0.001, 'start': 0, 'end': 10000000}
        self.assertEqual(apply_filters(df, min_blocks=5, min_proportion=0.01, block_size=1000000), 'low_proportion_me')

        # too short only
        df = {'mean_proportion_me': 0.5, 'start': 0, 'end': 2000000}
        self.assertEqual(apply_filters(df, min_blocks=5, min_proportion=0.01, block_size=1000000), 'min_block_length')

        # both
        df = {'mean_proportion_me': 0.001, 'start': 0, 'end': 2000000}
        self.assertEqual(apply_filters(df, min_blocks=5, min_proportion=0.01, block_size=1000000), 'low_proportion_me;min_block_length')

    def test_merge_contiguous_blocks(self):

        # empty df
        self.assertEqual(merge_contiguous_blocks(pd.DataFrame(), block_size=1000000, analysis='p_value', analysis2='proportion'), [])

        # single contiguous run on one chrom
        df = pd.DataFrame({
            'chrom': ['1', '1', '1'],
            'end': [1000000, 2000000, 3000000],
            'p_value': [0.1, 0.2, 0.3],
            'proportion': [0.5, 0.6, 0.7],
        })
        blocks = merge_contiguous_blocks(df, block_size=1000000, analysis='p_value', analysis2='proportion')
        self.assertEqual(len(blocks), 1)
        self.assertEqual(blocks[0][0], '1')       # chrom
        self.assertEqual(blocks[0][2], 3000000)   # last_end

        # gap within same chromosome -> splits into two blocks
        df = pd.DataFrame({
            'chrom': ['1', '1', '1'],
            'end': [1000000, 2000000, 5000000],  # gap between 2M and 5M
            'p_value': [0.1, 0.2, 0.3],
            'proportion': [0.5, 0.6, 0.7],
        })
        blocks = merge_contiguous_blocks(df, block_size=1000000, analysis='p_value', analysis2='proportion')
        self.assertEqual(len(blocks), 2)

        # chromosome change -> splits into two blocks
        df = pd.DataFrame({
            'chrom': ['1', '1', '2'],
            'end': [1000000, 2000000, 1000000],
            'p_value': [0.1, 0.2, 0.3],
            'proportion': [0.5, 0.6, 0.7],
        })
        blocks = merge_contiguous_blocks(df, block_size=1000000, analysis='p_value', analysis2='proportion')
        self.assertEqual(len(blocks), 2)
        self.assertEqual(blocks[0][0], '1')
        self.assertEqual(blocks[1][0], '2')

    from unittest.mock import MagicMock

    def test_get_variant_match_flags(self):

        variant = MagicMock()
        variant.matches_paternal_uniparental_ambiguous.return_value = True
        variant.matches_maternal_uniparental_ambiguous.return_value = False
        variant.matches_paternal_uniparental_isodisomy.return_value = False
        variant.matches_maternal_uniparental_isodisomy.return_value = False
        variant.alleles_identical_to_dad.return_value = True
        variant.alleles_identical_to_mum.return_value = False
        variant.is_hom_alt.return_value = True
        variant.is_biparental_inheritance.return_value = False

        flags = get_variant_match_flags(variant, proband_id='PROBAND1', min_gq=20, min_dp=20)
        self.assertTrue(flags['matches_paternal_uniparental_ambiguous_count'])
        self.assertTrue(flags['is_homozygous_count'])
        self.assertFalse(flags['is_biparental_count'])
        variant.is_hom_alt.assert_called_once_with('PROBAND1')

    def test_check_affected_and_sex(self):
        """
        Basic tests to ensure that this is exiting on unaffected and no sex samples
        """

        #Affected, sex set
        affected = '2'
        sex = '1'
        unaffected_or_no_sex = check_affected_and_sex(affected, sex)
        self.assertFalse(unaffected_or_no_sex)

        #Unaffected, sex set
        affected = '1'
        sex = '1'
        unaffected_or_no_sex = check_affected_and_sex(affected, sex)
        self.assertTrue(unaffected_or_no_sex)

        #Affected, no sex
        affected = '2'
        sex = '0'
        unaffected_or_no_sex = check_affected_and_sex(affected, sex)
        self.assertTrue(unaffected_or_no_sex)

        #Unaffected, no sex
        affected = '1'
        sex = '0'
        unaffected_or_no_sex = check_affected_and_sex(affected, sex)
        self.assertTrue(unaffected_or_no_sex)

if __name__ == "__main__":
    unittest.main()