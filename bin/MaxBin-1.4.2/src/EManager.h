#ifndef __EMANAGER_H__
#define __EMANAGER_H__

#include "global_inc.h"
#include "fastaReader.h"
#include "EucDist.h"
#include "kmerMap.h"
#include "Profiler.h"
#include "AbundanceLoader.h"
#include "logger.h"
#include <time.h>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/poisson.hpp>
using boost::math::normal;
using boost::math::poisson;

class EManager
{
	public:
		EManager(char *input_fs, char *input_abundance, char *input_output);
		~EManager();
		void run(char *seedfile);
		void setVerbose(bool is_verbose);
		void setMaxEM(int maxem);
		void setMinLength(unsigned int min_length);

	private:
		// Variables
		char *fastafile;
		char *outputfile;
		int kmer_len, kmer_len2;
		kmerMap *kmap, *kmap2;
		fastaReader *seq;
		AbundanceLoader *ab_loader;
		int max_EM;
		normal *intranormaldistr, *internormaldistr;
		normal *intranormaldistr2, *internormaldistr2;
		poisson **poissondistr;
		logger *log;
		time_t start_t, end_t;

		unsigned int min_seq_length;
		int seqnum;
		int seed_num;
		Profiler **seed_profile;
		Profiler **seq_profile;
		Profiler **seed_profile2;
		Profiler **seq_profile2;
		bool *is_profile_N;
		bool *is_profile2_N;
		double *seed_abundance;
		char **seed_header;
		double **seq_prob;
		int *seq_bin;
		int *seed_count;
		char **bin_name;
		bool *is_estimated;

		char str[1024];
		double MIN_PROB_THRESHOLD;
		int MAX_DIST_ESTIMATE;
		int FASTA_LINE;
		double VERY_SMALL_DOUBLE;
		int STABLE_BIN_COUNT;

		// Functions
		void init();
		void estimate_normaldistr();
		void init_EM();
		void run_EM(int run_time);
		bool classify(double min_prob, unsigned int min_seqlen);
		void filter_seed();
		void write_result();
		void write_fasta(char *seq, fstream *fs);
		double get_probability(double distance, double curr_abund, poisson *poisson_distr);
		double get_probability(double distance1, double distance2, double curr_abund, poisson *poisson_distr);
		double get_prob_dist(double distance);
		double get_prob_dist2(double distance);
		double get_prob_abund(double curr_abund, poisson *poisson_distr);
		void logtime_start();
		void logtime_end();
};


#endif
