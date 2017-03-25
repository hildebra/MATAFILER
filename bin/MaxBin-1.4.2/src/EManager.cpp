#include "EManager.h"
#include "SpearmanDist.h"
#include "ManhattanDist.h"
#include "EucDist.h"
#include <fstream>
#include <math.h>
#include "boost/random.hpp"
#include "boost/generator_iterator.hpp"

//#define USE_TWO_DIST

EManager::EManager(char *input_fs, char *input_abundance, char *input_output)
{
	logtime_start();
	fastafile = input_fs;
	seq = new fastaReader(input_fs);
	seqnum = seq->getNum();
	ab_loader = new AbundanceLoader(input_abundance);
	outputfile = input_output;
	init();
}

EManager::~EManager()
{
	int i, j;
	if (seed_profile != NULL)
	{
		for (i = 0; i < seed_num; i++)
		{
			delete(seed_profile[i]);
#ifdef USE_TWO_DIST
			delete(seed_profile2[i]);
#endif
			//free(seed_header[i]);
			if (seed_count[i] > 0)
			{
				free(bin_name[i]);
			}
		}
		free(seed_profile);
#ifdef USE_TWO_DIST
		free(seed_profile2);
#endif
		free(seed_header);
		free(seed_count);
		free(bin_name);
		j = seq->getNum();
		for (i = 0; i < j; i++)
		{
			delete(seq_profile[i]);
#ifdef USE_TWO_DIST
			delete(seq_profile2[i]);
#endif
			free(seq_prob[i]);
		}
		free(seq_profile);
		free(is_profile_N);
#ifdef USE_TWO_DIST
		free(seq_profile2);
		free(is_profile2_N);
#endif
		free(seq_prob);
		free(seed_abundance);
		free(seq_bin);
		free(is_estimated);
	}

	delete(kmap);
#ifdef USE_TWO_DIST
	delete(kmap2);
#endif
	if (intranormaldistr != NULL)
	{
		delete(intranormaldistr);
	}
	if (internormaldistr != NULL)
	{
		delete(internormaldistr);
	}
#ifdef USE_TWO_DIST
	if (intranormaldistr2 != NULL)
	{
		delete(intranormaldistr2);
	}
	if (internormaldistr2 != NULL)
	{
		delete(internormaldistr2);
	}
#endif
	delete(log);
	delete(seq);
	delete(ab_loader);
}

void EManager::init()
{
	seed_profile = NULL;
	seq_profile = NULL;
	is_profile_N = NULL;
#ifdef USE_TWO_DIST
	seed_profile2 = NULL;
	seq_profile2 = NULL;
	is_profile2_N = NULL
#endif
	seed_abundance = NULL;
	seq_prob = NULL;
	seq_bin = NULL;
	seed_count = NULL;
	bin_name = NULL;
	is_estimated = NULL;
	intranormaldistr = NULL;
	internormaldistr = NULL;
	intranormaldistr2 = NULL;
	internormaldistr2 = NULL;
	poissondistr = NULL;
	kmer_len = 4;
#ifdef USE_TWO_DIST
	kmer_len2 = 4;
#endif
	kmap = new kmerMap(kmer_len, true);
#ifdef USE_TWO_DIST
	kmap2 = new kmerMap(kmer_len2, true);
#endif
	min_seq_length = 1000;
	max_EM = 50;
	VERY_SMALL_DOUBLE = 1e-100;
	MIN_PROB_THRESHOLD = 0.8;
	MAX_DIST_ESTIMATE = 100000;
	FASTA_LINE = 70;
	STABLE_BIN_COUNT = 5;

	// Setup log file
	sprintf(str, "%s.log", outputfile);
	log = new logger(str);
	log->setVerbose(false);
}

void EManager::setVerbose(bool is_verbose)
{
	log->setVerbose(is_verbose);
}

void EManager::setMaxEM(int maxem)
{
	max_EM = maxem;
}

void EManager::setMinLength(unsigned int min_length)
{
	min_seq_length = min_length;
	sprintf(str, "Minimum contig length set to %d.\n", min_seq_length);
	log->writelog(str, true);
}

void EManager::run(char *seedfile)
{
	// Open seed file and look for the sequences used as seeds. One sequence per line
	fstream *fs = new fstream(seedfile, ios::in);

	sprintf(str, "Reading seed list...\n");
	log->writelog(str, true);

	seq->resetMark();
	seed_num = 0;
	while (!fs->eof())
	{
		memset(str, '\0', 1024);
		fs->getline(str, 1023);
		while (str[strlen(str) - 1] == '\n' || str[strlen(str) - 1] == '\r')
		{
			str[strlen(str) - 1] = '\0';
		}
		if (str[0] == '\0')
		{
			continue;
		}
		seq->setMark(str);
		seed_num++;

		log->writelog("\t", false);
		log->writelog(str, false);
		log->writelog("\n", false);

		if (fs->eof())
		{
			break;
		}
		else if (fs->fail())
		{
			fs->clear();
		}
	}
	fs->close();
	delete(fs);

	if (seed_num > 1)
	{
		init_EM();
		//sprintf(str, "Test run...\n");
		//log->writelog(str, true);
		//run_EM(10);
		// Only classify sequences >= threshold. For eliminating bins without sequences
		//classify(MIN_PROB_THRESHOLD, min_seq_length);
		//filter_seed();
		run_EM(max_EM);
		classify(MIN_PROB_THRESHOLD, min_seq_length);
		write_result();
		logtime_end();
	}
}

void EManager::init_EM()
{
	int i, j;
	double d;
	// Initialize EM data structures

	sprintf(str, "Looking for seeds in sequences.\n");
	log->writelog(str, true);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		if (seq->isMark(i) == true)
		{
			j++;
		}
	}
	if (j != seed_num)
	{
		sprintf(str, "Cannot find all seeds in contigs. Found %i seeds instead of %i in seed file.\nPlease check if the seed files.\n", j, seed_num);
		log->writelog(str, true);
		exit(-1);
	}
	seed_num = j;

	seq_profile = (Profiler**)malloc(sizeof(Profiler*) * seqnum);
	memset(seq_profile, '\0', sizeof(Profiler*) * seqnum);
	seed_profile = (Profiler**)malloc(sizeof(Profiler*) * seed_num);
	memset(seed_profile, '\0', sizeof(Profiler*) * seed_num);
	is_profile_N = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_profile_N, '\0', sizeof(bool) * seqnum);
#ifdef USE_TWO_DIST
	seq_profile2 = (Profiler**)malloc(sizeof(Profiler*) * seqnum);
	memset(seq_profile2, '\0', sizeof(Profiler*) * seqnum);
	seed_profile2 = (Profiler**)malloc(sizeof(Profiler*) * seed_num);
	memset(seed_profile2, '\0', sizeof(Profiler*) * seed_num);
	is_profile2_N = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_profile2_N, '\0', sizeof(bool) * seqnum);
#endif
	seed_abundance = (double*)malloc(sizeof(double) * seed_num);
	memset(seed_abundance, '\0', sizeof(double) * seed_num);
	seed_header = (char**)malloc(sizeof(char*) * seed_num);
	memset(seed_header, '\0', sizeof(char*) * seed_num);
	seed_count = (int*)malloc(sizeof(int) * seed_num);
	memset(seed_count, '\0', sizeof(int) * seed_num);
	seq_bin = (int*)malloc(sizeof(int) * seqnum);
	memset(seq_bin, '\0', sizeof(int) * seqnum);
	is_estimated = (bool*)malloc(sizeof(bool) * seqnum);
	memset(is_estimated, '\0', sizeof(bool) * seqnum);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		seq_profile[i] = new Profiler(kmer_len, seq->getSeqByNum(i), kmap);
		if (seq_profile[i]->getPercentN() == 1)
		{
			is_profile_N[i] = true;
		}
#ifdef USE_TWO_DIST
		seq_profile2[i] = new Profiler(kmer_len2, seq->getSeqByNum(i), kmap2);
		if (seq_profile2[i]->getPercentN() == 1)
		{
			is_profile2_N[i] = true;
		}
#endif
		d = ab_loader->getAbundance(seq->getHeaderByNum(i));
		if (seq->isMark(i) == true && d > 0)
		{
			seed_profile[j] = new Profiler(kmer_len, seq->getSeqByNum(i), kmap);
#ifdef USE_TWO_DIST
			seed_profile2[j] = new Profiler(kmer_len2, seq->getSeqByNum(i), kmap2);
#endif
			//seed_abundance[j] = ab_loader->getAbundance(seq->getHeaderByNum(i));
			seed_abundance[j] = d;
			seed_header[j] = seq->getHeaderByNum(i);
			sprintf(str, "\t%s [%f]\n", seed_header[j], d);
			log->writelog(str, true);
			j++;
		}
		else
		{
			seq->isMark(i) == false;
		}
	}
	estimate_normaldistr();

	if (j != seed_num)
	{
		seed_num = j;
	}
	sprintf(str, "Get %d seeds.\n", seed_num);
	log->writelog(str, true);

	/*
	double d;
	EucDist edist(kmer_len);
	for (i = 0; i < seed_num; i++)
	{
		for (j = i + 1; j < seed_num; j++)
		{
			d = edist.getDist(seed_profile[i]->getProfile(), seed_profile[j]->getProfile());
			if ((1 - cdf(*intranormaldistr, d)) >= 0.97)
			{
				printf("%s and %s are quite close at %f probability.\n", seed_header[i], seed_header[j], 1 - cdf(*intranormaldistr, d));
			}
		}
	}
	exit(-1);
	*/

	seq_prob = (double**)malloc(sizeof(double*) * seqnum);
	for (i = 0; i < seqnum; i++)
	{
		seq_prob[i] = (double*)malloc(sizeof(double) * seed_num);
		memset(seq_prob[i], '\0', sizeof(double) * seed_num);
	}
}

void EManager::estimate_normaldistr()
{
	double mean, std, mean2, std2;
	double mean_2, std_2, mean2_2, std2_2;
	// Pre-defined mean and standard deviation found from simulation test on 3181 IMG bacterial and archaeal genomes

	// Euclidean distance
	// 4-mer
	mean = 0;
	std = 0.01037897 / 2;
	//mean = 0.01505219 / 2;
	//std = 0.01037897 / 2;
	//mean = 0;
	//std = 0.01037897;
	mean2 = 0.0676654;
	std2 = 0.03419337;
	// 5-mer
	mean_2 = 0.01228021 / 2;
	std_2 = 0.008654644 / 2;
	//mean_2 = 0;
	//std_2 = 0.008654644;
	mean2_2 = 0.04346744;
	std2_2 = 0.02088791;
	// 6-mer
	//mean = 0.01032181;
	//std = 0.007806874;
	//mean = 0.02672312;
	//std = 0.01257017;

	// Spearman footrule distance
	// 4-mer
	//mean = 0.05849619;
	//std = 0.03473218;
	//mean2 = 0.2523695;
	//std2 = 0.113352;
	mean_2 = 0.05849619 / 2;
	std_2 = 0.03473218 / 2;
	mean2_2 = 0.2523695;
	std2_2 = 0.113352;

	// Manhattan distance
	// 4-mer
	//mean = 0.06337266;
	//std = 0.03575501;
	//mean2 = 0.2886038;
	//std2 = 0.1432431;

	intranormaldistr = new normal(mean, std);
	internormaldistr = new normal(mean2, std2);
#ifdef USE_TWO_DIST
	intranormaldistr2 = new normal(mean_2, std_2);
	internormaldistr2 = new normal(mean2_2, std2_2);
#endif

	// Estimate mean and variance for the distances in the input dataset.
	/*
	int i, n1, n2;
	double *dist;
	double d, mean, std;
	EucDist edist(kmer_len);
	//SpearmanDist edist(kmer_len);
	edist.setNormalization(true);

	// Boost random number generator
	boost::mt19937 rng;
	boost::uniform_int<> drawer(0, seq->getNum() - 1);
	boost::variate_generator<boost::mt19937, boost::uniform_int<> > dice(rng, drawer);
	sprintf(str, "Estimating the parameters for distance probability function.\nSampling %d pairs of sequences...\n", MAX_DIST_ESTIMATE);
	log->writelog(str, true);

	dist = (double*)malloc(sizeof(double) * MAX_DIST_ESTIMATE);
	memset(dist, '\0', sizeof(double) * MAX_DIST_ESTIMATE);

	mean = 0;
	for (i = 0; i < MAX_DIST_ESTIMATE; i++)
	{
		n1 = 0;
		n2 = 0;
		while (n1 == n2)
		{
			n1 = dice();
			n2 = dice();
		}
		d = edist.getDist(seq_profile[n1]->getProfile(), seq_profile[n2]->getProfile());
		dist[i] = d;
		mean = mean + d;
	}

	mean = mean / MAX_DIST_ESTIMATE; // Get mean value
	std = 0;
	for (i = 0; i < MAX_DIST_ESTIMATE; i++)
	{
		std = std + pow(mean - dist[i], 2);
	}
	std = sqrt(std / (MAX_DIST_ESTIMATE - 1));

	sprintf(str, "Mean = %f, Standard deviation = %f\n", mean, std);
	log->writelog(str, false);

	free(dist);

	// Create normal distribution
	intranormaldistr = new normal(mean, std);
	*/
}

void EManager::run_EM(int run_time)
{
	int run, i, j, k, nan, diff_count, tempbin, stable_count;
	double sum, d, d2, **seed_prob, max;
	double *dist_prob, *dist_prob2, *abund_prob;
	EucDist edist(kmer_len);
	edist.setNormalization(true);
#ifdef USE_TWO_DIST
	//EucDist edist2(kmer_len2);
	//SpearmanDist edist(kmer_len);
	SpearmanDist edist2(kmer_len2);
	//ManhattanDist edist(kmer_len);
	edist2.setNormalization(true);
#endif

	dist_prob = (double*)malloc(sizeof(double) * seed_num);
	memset(dist_prob, '\0', sizeof(double) * seed_num);
#ifdef USE_TWO_DIST
	dist_prob2 = (double*)malloc(sizeof(double) * seed_num);
	memset(dist_prob2, '\0', sizeof(double) * seed_num);
#endif
	abund_prob = (double*)malloc(sizeof(double) * seed_num);
	memset(abund_prob, '\0', sizeof(double) * seed_num);

	sprintf(str, "\nStart EM process.\n");
	log->writelog(str, true);

	// Execute EM algorithm
	stable_count = 0;
	for (run = 0; run < run_time; run++)
	{
		sprintf(str, "Iteration %d\n", run + 1);
		log->writelog(str, true);

		// Expectation
		// For each run, calculate the distance between the sequences and seeds and calculate probability
		poissondistr = (poisson**)malloc(sizeof(poisson*) * seed_num);
		memset(poissondistr, '\0', sizeof(poisson*) * seed_num);
		for (i = 0; i < seed_num; i++)
		{
			sprintf(str, "Abundance info [%s] [%d]: %f\n", seed_header[i], i + 1, seed_abundance[i]);
			log->writelog(str, false);
			if (seed_abundance[i] > 0)
			{
				poissondistr[i] = new poisson(seed_abundance[i]);
			}
			else
			{
				poissondistr[i] = NULL;
			}
		}
		nan = 0;
		diff_count = 0;
		for (i = 0; i < seqnum; i++)
		{
			if (seq->getSeqLenByNum(i) >= min_seq_length && is_profile_N[i] == false)
			{
				// Test of separating abundance probability and tetramer probability
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist(seq_profile[i]->getProfile(), seed_profile[j]->getProfile());
					dist_prob[j] = get_prob_dist(d);
					sum = sum + dist_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob[j] = dist_prob[j] / sum;
				}
#ifdef USE_TWO_DIST
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist2.getDist(seq_profile2[i]->getProfile(), seed_profile2[j]->getProfile());
					dist_prob2[j] = get_prob_dist2(d);
					sum = sum + dist_prob2[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob2[j] = dist_prob2[j] / sum;
				}
#endif
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					abund_prob[j] = get_prob_abund(ab_loader->getAbundance(seq->getHeaderByNum(i)), poissondistr[j]);
					sum = sum + abund_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					abund_prob[j] = abund_prob[j] / sum;
				}
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
#ifdef USE_TWO_DIST
					seq_prob[i][j] = dist_prob[j] * dist_prob2[j] * abund_prob[j];
#else
					seq_prob[i][j] = dist_prob[j] * abund_prob[j];
#endif
					sum = sum + seq_prob[i][j];
				}

				/*
				if (sum < VERY_SMALL_DOUBLE)
				{
					sum = 1;
				}
				*/
				k = 0;
				max = 0;
				tempbin = -1;
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				log->writelog(str, false);
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					if (max < seq_prob[i][j])
					{
						max = seq_prob[i][j];
						tempbin = j;
					}
					// Special handling of NAN (caused when sum is too small...)
					//sprintf(str, "\t(%d)%f [%f * %f * %f]", j + 1, seq_prob[i][j], dist_prob[j], dist_prob2[j], abund_prob[j]);
					sprintf(str, "\t(%d)%f", j + 1, seq_prob[i][j]);
					log->writelog(str, false);
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
						k = 1;
					}
				}
				if (k == 0 && seq_bin[i] != tempbin)
				{
					diff_count++;
					seq_bin[i] = tempbin;
				}
				log->writelog("\n", false);
				if (k == 1)
				{
					nan++;
				}
				is_estimated[i] = true;
			}
		}
		sprintf(str, "NAN: %d\n", nan);
		log->writelog(str, false);

		if (diff_count == 0)
		{
			stable_count++;
		}
		else
		{
			stable_count = 0;
		}

		// Maximization
		for (i = 0; i < seed_num; i++)
		{
			seed_abundance[i] = 0;
			seed_profile[i]->reset();
#ifdef USE_TWO_DIST
			seed_profile2[i]->reset();
#endif
			d = 0;
			for (j = 0; j < seqnum; j++)
			{
				if (seq->getSeqLenByNum(j) >= min_seq_length && seq_prob[j][i] == seq_prob[j][i] && is_profile_N[j] == false)
				{
					// Update seed abundances
					seed_abundance[i] = seed_abundance[i] + (ab_loader->getAbundance(seq->getHeaderByNum(j)) * seq->getSeqLenByNum(j) * seq_prob[j][i]);
					d = d + seq->getSeqLenByNum(j) * seq_prob[j][i];
					// Update seed profiles
					seed_profile[i]->addProfile(seq_profile[j], seq->getSeqLenByNum(j) * seq_prob[j][i]);
#ifdef USE_TWO_DIST
					seed_profile2[i]->addProfile(seq_profile2[j], seq->getSeqLenByNum(j) * seq_prob[j][i]);
#endif
				}
			}
			if (d == 0)
			{
				seed_abundance[i] = -1;
			}
			else
			{
				seed_profile[i]->calcProfile();
#ifdef USE_TWO_DIST
				seed_profile2[i]->calcProfile();
#endif
				seed_abundance[i] = seed_abundance[i] / d;
			}
		}
		/*
		for (i = 0; i < seed_num; i++)
		{
			if (seed_abundance[i] == -1)
			{
printf("Before %d: [i - 1] %f [i] %f [i + 1] %f\n", i + 1, seed_abundance[i - 1], seed_abundance[i], seed_abundance[i + 1]);
				delete(seed_profile[i]);
				delete(seed_profile2[i]);
				memcpy(&(seed_profile[i]), &(seed_profile[i + 1]), sizeof(Profiler*) * (seed_num - i - 1));
				memcpy(&(seed_profile2[i]), &(seed_profile2[i + 1]), sizeof(Profiler*) * (seed_num - i - 1));
				memcpy(&(seed_abundance[i]), &(seed_abundance[i + 1]), sizeof(double) * (seed_num - i - 1));
				memcpy(&(seed_header[i]), &(seed_header[i + 1]), sizeof(char*) * (seed_num - i - 1));
				seed_num--;
				i--;
				sprintf(str, "All contigs assigned to bin %d with 0 probability. Remove the bin.\n", i + 1);
				log->writelog(str, false);
printf("After %d: [i - 1] %f [i] %f [i + 1] %f\n", i + 1, seed_abundance[i - 1], seed_abundance[i], seed_abundance[i + 1]);
			}
		}
		*/

		for (i = 0; i < seed_num; i++)
		{
			if (poissondistr[i] != NULL)
			{
				delete(poissondistr[i]);
			}
		}
		free(poissondistr);

		if (stable_count >= STABLE_BIN_COUNT)
		{
			break;
		}
	}
	sprintf(str, "\nEM finishes successfully.\n");
	log->writelog(str, true);

	sprintf(str, "%s.dist", outputfile);
	fstream *dist = new fstream(str, ios::out);
	poissondistr = (poisson**)malloc(sizeof(poisson*) * seed_num);
	memset(poissondistr, '\0', sizeof(poisson*) * seed_num);
	seed_prob = (double**)malloc(sizeof(double*) * seed_num);
	memset(seed_prob, '\0', sizeof(double*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		if (seed_abundance[i] > 0)
		{
			poissondistr[i] = new poisson(seed_abundance[i]);
		}
		else
		{
			poissondistr[i] = NULL;
		}
		seed_prob[i] = (double*)malloc(sizeof(double) * seed_num);
		memset(seed_prob[i], '\0', sizeof(double) * seed_num);
	}
	for (i = 0; i < seed_num; i++)
	{
		sprintf(str, "\tBin%3d", i + 1);
		dist->write(str, strlen(str));
	}
	dist->write("\n", 1);
	for (i = 0; i < seed_num; i++)
	{
		sum = 0;
		for (j = 0; j < seed_num; j++)
		{
			d = edist.getDist(seed_profile[i]->getProfile(), seed_profile[j]->getProfile());
#ifdef USE_TWO_DIST
			d2 = edist2.getDist(seed_profile2[i]->getProfile(), seed_profile2[j]->getProfile());
			seed_prob[i][j] = get_probability(d, d2, seed_abundance[j], poissondistr[i]);
#else
			seed_prob[i][j] = get_probability(d, seed_abundance[j], poissondistr[i]);
#endif
			sum = sum + seed_prob[i][j];
		}
		sprintf(str, "Bin%3d", i + 1);
		dist->write(str, strlen(str));
		for (j = 0; j < seed_num; j++)
		{
			seed_prob[i][j] = seed_prob[i][j] / sum;
			// Special handling of NAN (caused when sum is too small...)
			if (seed_prob[i][j] != seed_prob[i][j])
			{
				seed_prob[i][j] = 0;
			}
			sprintf(str, "\t%f", seed_prob[i][j]);
			dist->write(str, strlen(str));
		}
		dist->write("\n", 1);
	}
	dist->close();
	delete(dist);
	for (i = 0; i < seed_num; i++)
	{
		if (poissondistr[i] != NULL)
		{
			delete(poissondistr[i]);
		}
		free(seed_prob[i]);
	}
	free(seed_prob);
	free(poissondistr);
	free(dist_prob);
#ifdef USE_TWO_DIST
	free(dist_prob2);
#endif
	free(abund_prob);
}

bool EManager::classify(double min_prob, unsigned int min_seqlen)
{
	int i, j;
	double max, sum, d;
	double *dist_prob, *dist_prob2, *abund_prob;
	bool ret;
	poisson **poissondistr;

	sprintf(str, "%s.prob_dist", outputfile);
	fstream *distf = new fstream(str, ios::out);

	EucDist edist(kmer_len);
	edist.setNormalization(true);
#ifdef USE_TWO_DIST
	//EucDist edist2(kmer_len2);
	//SpearmanDist edist(kmer_len);
	SpearmanDist edist2(kmer_len2);
	//ManhattanDist edist(kmer_len);
	edist2.setNormalization(true);
#endif

	dist_prob = (double*)malloc(sizeof(double) * seed_num);
	memset(dist_prob, '\0', sizeof(double) * seed_num);
#ifdef USE_TWO_DIST
	dist_prob2 = (double*)malloc(sizeof(double) * seed_num);
	memset(dist_prob2, '\0', sizeof(double) * seed_num);
#endif
	abund_prob = (double*)malloc(sizeof(double) * seed_num);
	memset(abund_prob, '\0', sizeof(double) * seed_num);

	sprintf(str, "\nClassifying sequences based on the EM result.\nMinimum probability for binning: %0.2f\n", MIN_PROB_THRESHOLD);
	log->writelog(str, true);

	poissondistr = (poisson**)malloc(sizeof(poisson*) * seed_num);
	memset(poissondistr, '\0', sizeof(poisson*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		if (seed_abundance[i] > 0)
		{
			poissondistr[i] = new poisson(seed_abundance[i]);
		}
		else
		{
			poissondistr[i] = NULL;
		}
	}

	for (i = 0; i < seqnum; i++)
	{
		if (seq->getSeqLenByNum(i) >= min_seqlen && is_profile_N[i] == false)
		{
			if (is_estimated[i] == false)
			{
				// Test of separating abundance probability and tetramer probability
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist(seq_profile[i]->getProfile(), seed_profile[j]->getProfile());
					dist_prob[j] = get_prob_dist(d);
					sum = sum + dist_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob[j] = dist_prob[j] / sum;
				}
#ifdef USE_TWO_DIST
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist2.getDist(seq_profile2[i]->getProfile(), seed_profile2[j]->getProfile());
					dist_prob2[j] = get_prob_dist2(d);
					sum = sum + dist_prob2[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					dist_prob2[j] = dist_prob2[j] / sum;
				}
#endif
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					abund_prob[j] = get_prob_abund(ab_loader->getAbundance(seq->getHeaderByNum(i)), poissondistr[j]);
					sum = sum + abund_prob[j];
				}
				for (j = 0; j < seed_num; j++)
				{
					abund_prob[j] = abund_prob[j] / sum;
				}
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
#ifdef USE_TWO_DIST
					seq_prob[i][j] = dist_prob[j] * dist_prob2[j] * abund_prob[j];
#else
					seq_prob[i][j] = dist_prob[j] * abund_prob[j];
#endif
					sum = sum + seq_prob[i][j];
				}
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				log->writelog(str, false);
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%f", j + 1, seq_prob[i][j]);
					log->writelog(str, false);
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				log->writelog("\n", false);
			}

			max = 0;
			for (j = 0; j < seed_num; j++)
			{
				if (seq_prob[i][j] > max)
				{
					max = seq_prob[i][j];
					seq_bin[i] = j;
				}
			}
			if (max <= min_prob)
			{
				seq_bin[i] = -1;
			}
			else
			{
				seed_count[seq_bin[i]]++;
			}
			if (seq_bin[i] != -1)
			{
				sprintf(str, "seq [%s] assigned to bin [%d]\n", seq->getHeaderByNum(i), seq_bin[i]);
				log->writelog(str, false);
			}

				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					d = edist.getDist(seq_profile[i]->getProfile(), seed_profile[j]->getProfile());
					seq_prob[i][j] = get_prob_dist(d);
					sum = sum + seq_prob[i][j];
				}
				sprintf(str, "[%s]", seq->getHeaderByNum(i));
				distf->write(str, strlen(str));
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%f", j + 1, seq_prob[i][j]);
					distf->write(str, strlen(str));
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				distf->write("\n", 1);
				sum = 0;
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = get_prob_abund(ab_loader->getAbundance(seq->getHeaderByNum(i)), poissondistr[j]);
					sum = sum + seq_prob[i][j];
				}
				for (j = 0; j < seed_num; j++)
				{
					seq_prob[i][j] = seq_prob[i][j] / sum;
					// Special handling of NAN (caused when sum is too small...)
					sprintf(str, "\t(%d)%f", j + 1, seq_prob[i][j]);
					distf->write(str, strlen(str));
					if (seq_prob[i][j] != seq_prob[i][j])
					{
						seq_prob[i][j] = 0;
					}
				}
				distf->write("\n", 1);
		}
		else
		{
			seq_bin[i] = -1;
		}
	}

	for (i = 0; i < seed_num; i++)
	{
		if (poissondistr[i] != NULL)
		{
			delete(poissondistr[i]);
		}
	}
	free(poissondistr);

	distf->close();
	delete(distf);
	free(dist_prob);
#ifdef USE_TWO_DIST
	free(dist_prob2);
#endif
	free(abund_prob);

	ret = false;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] == 0)
		{
			ret = true;
			break;
		}
	}
	return(ret);
}

void EManager::filter_seed()
{
	// seed_profile
	// seed_abundance
	// seed_header
	// seed_num
	int i, j;
	int temp_num;
	Profiler **temp_profile;
	double *temp_abundance;
	char **temp_header;
	int *temp_count;

	// Get the number of bins that are not empty
	temp_num = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			temp_num++;
		}
	}
	sprintf(str, "Filter out %d bins with no results during test run.\n", seed_num - temp_num);
	log->writelog(str, true);

	temp_profile = (Profiler**)malloc(sizeof(Profiler*) * temp_num);
	memset(temp_profile, '\0', sizeof(Profiler*) * temp_num);
	temp_abundance = (double*)malloc(sizeof(double) * temp_num);
	memset(temp_abundance, '\0', sizeof(double) * temp_num);
	temp_header = (char**)malloc(sizeof(char*) * temp_num);
	memset(temp_header, '\0', sizeof(char*) * temp_num);
	temp_count = (int*)malloc(sizeof(int) * temp_num);
	memset(temp_count, '\0', sizeof(int) * temp_num);

	j = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			temp_profile[j] = seed_profile[i];
			temp_abundance[j] = seed_abundance[i];
			temp_header[j] = seed_header[i];
			temp_count[j] = seed_count[i];
			j++;
		}
		else
		{
			delete(seed_profile[i]);
		}
	}
	free(seed_profile);
	free(seed_abundance);
	free(seed_header);
	free(seed_count);
	seed_profile = temp_profile;
	seed_abundance = temp_abundance;
	seed_header = temp_header;
	seed_num = temp_num;
	seed_count = temp_count;

	for (i = 0; i < seqnum; i++)
	{
		free(seq_prob[i]);
	}
	free(seq_prob);

	seq_prob = (double**)malloc(sizeof(double*) * seqnum);
	for (i = 0; i < seqnum; i++)
	{
		seq_prob[i] = (double*)malloc(sizeof(double) * seed_num);
		memset(seq_prob[i], '\0', sizeof(double) * seed_num);
	}
}

void EManager::write_result()
{
	int i, j, k;
	fstream **fs, *un;

	// Name the bins -- ignore bins with no sequences.
	bin_name = (char**)malloc(sizeof(char*) * seed_num);
	memset(bin_name, '\0', sizeof(char*) * seed_num);
	j = 1;
	k = 0;
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			bin_name[i] = (char*)malloc(sizeof(char) * (strlen(outputfile) + 20));
			memset(bin_name[i], '\0', sizeof(char) * (strlen(outputfile) + 20));
			sprintf(bin_name[i], "%s.%03d.fasta", outputfile, j);
			j++;
		}
		else
		{
			k++;
		}
	}
	sprintf(str, "Ignoring %d bins without any sequences.\n", k);
	log->writelog(str, true);

	// Write summary report into report page
	sprintf(str, "%s.summary", outputfile);
	un = new fstream(str, ios::out);
	for(i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			sprintf(str, "Bin [%s] (%f)\n", bin_name[i],seed_abundance[i]);
			un->write(str, strlen(str));
			for (j = 0; j < seqnum; j++)
			{
				if (seq_bin[j] == i)
				{
					sprintf(str, "\t%s (%f)\n", seq->getHeaderByNum(j), ab_loader->getAbundance(seq->getHeaderByNum(j)));
					un->write(str, strlen(str));
				}
			}
			un->write("\n", 1);
		}
	}
	sprintf(str, "\nBins without any sequences:\n");
	un->write(str, strlen(str));
	j = 1;
	for(i = 0; i < seed_num; i++)
	{
		if (seed_count[i] == 0)
		{
			sprintf(str, "%d: abundance %f\n", j, seed_abundance[i]);
			un->write(str, strlen(str));
			j++;
		}
	}
	un->close();
	delete(un);

	fs = (fstream**)malloc(sizeof(fstream*) * seed_num);
	memset(fs, '\0', sizeof(fstream*) * seed_num);
	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			fs[i] = new fstream(bin_name[i], ios::out);
		}
	}
	sprintf(str, "%s.noclass", outputfile);
	un = new fstream(str, ios::out);

	j = 0;
	for (i = 0; i < seqnum; i++)
	{
		if (seq_bin[i] != -1)
		{
			sprintf(str, ">%s\n", seq->getHeaderByNum(i));
			fs[seq_bin[i]]->write(str, strlen(str));
			write_fasta(seq->getSeqByNum(i), fs[seq_bin[i]]);
		}
		else
		{
			j++;
			sprintf(str, ">%s\n", seq->getHeaderByNum(i));
			un->write(str, strlen(str));
			write_fasta(seq->getSeqByNum(i), un);
		}
	}

	for (i = 0; i < seed_num; i++)
	{
		if (seed_count[i] > 0)
		{
			fs[i]->close();
			delete(fs[i]);
		}
	}
	free(fs);
	un->close();
	delete(un);

	sprintf(str, "Number of unclassified sequences: %d (%2.2f%%)\n", j, ((float)j / (float)(seq->getNum())) * 100);
	log->writelog(str, true);
}

void EManager::write_fasta(char *seq, fstream *fs)
{
	int len, i;
	char *c = seq;
	len = strlen(seq);
	i = 0;
	while(i < len)
	{
		if (i + FASTA_LINE >= len)
		{
			i = len;
			fs->write(c, strlen(c));
			fs->write("\n", 1);
		}
		else
		{
			i = i + FASTA_LINE;
			fs->write(c, FASTA_LINE);
			fs->write("\n", 1);
			c = c + FASTA_LINE;
		}
	}
}

double EManager::get_probability(double distance, double curr_abund, poisson *poisson_distr)
{
	double ret = get_prob_dist(distance) * get_prob_abund(curr_abund, poisson_distr);
	return ret;
}

double EManager::get_probability(double distance1, double distance2, double curr_abund, poisson *poisson_distr)
{
	double ret = get_prob_dist(distance1) * get_prob_dist(distance2) * get_prob_abund(curr_abund, poisson_distr);
	return ret;
}

double EManager::get_prob_dist(double distance)
{
	//double ret = (1 - cdf(*intranormaldistr, distance));
	double d_intra = pdf(*intranormaldistr, distance);
	double d_inter = pdf(*internormaldistr, distance);
	double ret = d_intra / (d_inter + d_intra);
	//double ret = pow(d_intra / (d_inter + d_intra), 2);
	return ret;
}

/* Using hardwired distance distribution -- not performing well
double EManager::get_prob_dist(double distance)
{
	// Initialize hardwired distribution
	double dist_prob[650] = {0.0457,0.0381,0.0422,0.0402,0.0438,0.0391,0.0411,0.0372,0.0372,0.0331,0.0339,0.0334,0.0298,0.0290,0.0265,0.0262,0.0249,0.0218,0.0215,0.0207,0.0189,0.0171,0.0186,0.0149,0.0160,0.0142,0.0123,0.0113,0.0106,0.0106,0.0106,0.0092,0.0078,0.0092,0.0083,0.0065,0.0067,0.0064,0.0065,0.0062,0.0058,0.0051,0.0050,0.0043,0.0045,0.0042,0.0035,0.0033,0.0036,0.0029,0.0029,0.0029,0.0028,0.0024,0.0026,0.0024,0.0023,0.0023,0.0017,0.0020,0.0018,0.0016,0.0018,0.0016,0.0012,0.0013,0.0015,0.0015,0.0011,0.0009,0.0011,0.0008,0.0010,0.0009,0.0010,0.0010,0.0007,0.0007,0.0008,0.0008,0.0006,0.0007,0.0007,0.0006,0.0006,0.0005,0.0006,0.0005,0.0005,0.0006,0.0004,0.0004,0.0005,0.0003,0.0005,0.0004,0.0004,0.0002,0.0003,0.0002,0.0002,0.0003,0.0002,0.0001,0.0003,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0001,0.0002,0.0002,0.0002,0.0001,0,0.0001,0.0002,0.0002,0.0002,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0.0001,0,0,0.0001,0.0001,0.0001,0.0001,0,0.0001,0,0,0,0,0,0,0.0001,0.0001,0.0001,0,0,0,0,0,0.0001,0,0,0,0,0,0.0001,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0001,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	double interval = 0.001;

	int i = distance / interval;
	if (i >= 650)
	{
		return 0;
	}
	else
	{
		return dist_prob[i], 2;
	}
}
*/

double EManager::get_prob_dist2(double distance)
{
	//double ret = (1 - cdf(*intranormaldistr2, distance));
	double d_intra = pdf(*intranormaldistr2, distance);
	double d_inter = pdf(*internormaldistr2, distance);
	double ret = d_intra / (d_inter + d_intra);
	return ret;
}

double EManager::get_prob_abund(double curr_abund, poisson *poisson_distr)
{
	double ret;
	if (poisson_distr != NULL && curr_abund > 0)
	{
		//ret = sqrt(pdf(*poisson_distr, curr_abund));
		ret = pdf(*poisson_distr, curr_abund);
	}
	else
	{
		ret = 0;
	}
	return ret;
}

void EManager::logtime_start()
{
	start_t = time(NULL);
}

void EManager::logtime_end()
{
	end_t = time(NULL);
	log->writetime(end_t - start_t);
}
