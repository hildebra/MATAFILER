#include "EManager.h"

int main(int argc, char *argv[])
{
	int i;
	char *inputfasta = NULL, *seed_file = NULL, *abund_file = NULL, *out = NULL;
	int maxem = -1;
	bool free_out = false;
	bool isverbose = false;
	unsigned int min_length = 0;

	// Handle parameters
	for (i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-fasta") == 0)
		{
			i++;
			inputfasta = argv[i];
		}
		else if (strcmp(argv[i], "-out") == 0)
		{
			i++;
			out = argv[i];
		}
		else if (strcmp(argv[i], "-seed") == 0)
		{
			i++;
			seed_file = argv[i];
		}
		else if (strcmp(argv[i], "-abund") == 0)
		{
			i++;
			abund_file = argv[i];
		}
		else if (strcmp(argv[i], "-max_run") == 0)
		{
			i++;
			maxem = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-min_contig_length") == 0)
		{
			i++;
			min_length = atoi(argv[i]);
		}
		else if (strcmp(argv[i], "-verbose") == 0)
		{
			isverbose = true;
		}
		else
		{
			printf("Unrecognized token: %s\n", argv[i]);
			printf("Usage: MaxBin -fasta (input fasta) -seed (seed file) -abund (abundance file) [-out (output file)]\n");
			exit(-1);
		}
	}

	if (inputfasta == NULL || seed_file == NULL || abund_file == NULL)
	{
		printf("Usage: MaxBin -fasta (input fasta) -seed (seed file) -abund (abundance file) [-out (output file)] [-max_run (max EM iteration; default 50)]\n");
		exit(-1);
	}
	if (out == NULL)
	{
		out = (char*)malloc(sizeof(char) * (strlen(inputfasta) + 5));
		memset(out, '\0', sizeof(char) * (strlen(inputfasta) + 5));
		sprintf(out, "%s.out", inputfasta);
		free_out = true;
	}

	EManager em(inputfasta, abund_file, out);
	em.setVerbose(isverbose);
	if (maxem != -1)
	{
		em.setMaxEM(maxem);
	}
	if (min_length > 0)
	{
		em.setMinLength(min_length);
	}
	em.run(seed_file);

	if (free_out == true)
	{
		free(out);
	}
}

