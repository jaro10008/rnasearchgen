#include<cstdio>
#include<string>
#include<cstdlib>
#include<ctime>
#include<cmath>
#include<algorithm>
#include<omp.h>
#include<vector>
#include<cctype>

extern "C" {
	#include<ViennaRNA/fold.h>
}

using namespace std;

static char nucToChar(int nuc){
	switch(nuc){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'U';
		default:
			return '?';
	}
}

static int charToNuc(char ch){
	switch(ch){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'U':
			return 3;

		default:
			return -1;
	}
}

static bool checkIfMatch(char ch, unsigned char temp){
	switch(ch){
		case 'A':
			return bool(temp & 1);
		case 'C':
			return bool(temp & 2);
		case 'G':
			return bool(temp & 4);
		case 'U':
			return bool(temp & 8);
		case '(':
			return bool(temp & 16);
		case ')':
			return bool(temp & 32);
		case '.':
			return bool(temp & 64);
		default:
			return false;
	}
}

static vector<unsigned char> temp;
static bool begi, endi;


struct sek{
	string seq;
	string str;
	float mfe;
	int result;

	sek(unsigned int ss = 0){
		if(ss){
			for(int i = 0; i < ss; ++i){
				int nuk = rand()%4;
				seq += nucToChar(nuk);
			}
		
			run_fold();

			result = rate();
		}
	}
	sek(const char* src){
		seq += src;

		run_fold();

		result = rate();
	}

	bool operator<(sek a) const{
		if(result == a.result && result == 0.0)
			return mfe < a.mfe;
		return result < a.result;
	}

	private:

	int rate(){
		unsigned int s = seq.size();
		unsigned int ts = temp.size();

		if(s < ts){
			return int(2 * ts - s);
		}
		else if(begi && endi && ts != s)
			return int(s);

		unsigned int minpos = 0;
		unsigned int maxpos = s - ts;

		if(begi){
			maxpos = 0;
		}
		if(endi){
			minpos = s - ts;
		}

		unsigned int minRes = ts;

		for(unsigned int i = minpos; i <= maxpos; ++i){
			unsigned int res = 0;
			for(unsigned int j = 0; j < ts; ++j){
				bool n1 = checkIfMatch(seq[i + j], temp[j]);
				bool n2 = checkIfMatch(str[i + j], temp[j]);

				if(!n1){
					++res;
				}
				if(!n2){
					++res;
				}
			}

			if(res < minRes){
				minRes = res;
			}
		}

		return int(minRes);
	}

	void run_fold(){
		char* temp = new char[seq.size() + 1];
		
		mfe = vrna_fold(seq.c_str(), temp);
		str = temp;

		delete []temp;
	}
};

static const int M = 1000;

static sek tab[M], tab2[M];
static sek* curr = tab;
static sek* curr2 = tab2;

static unsigned int randomNumber(unsigned int n){
	unsigned int r1 = rand();
	unsigned int r2 = rand();
	unsigned int r3 = rand();

	return ((r1 << 30) ^ (r2 << 15) ^ r3) % n;
}

static int b = 100;

static int choose(){/*
	int in = (int) sqrt(randomNumber(M * M));
	return (M - 1 - in);*/

	int a = 1;
	int p = randomNumber(M * (a * M + b));
	int in = (sqrt(b * b + 4 * a * p) - b) / (2*a);

	return  M - 1 - in;

}

static int mutation = 256;

static sek cross(const sek& s1, const sek& s2){
	string seq;

	int in = 0;
	bool p = (bool) rand() % 2;

	while(in < s1.seq.size() && in < s2.seq.size()){

		if(rand()%(mutation) == 0){
			seq += nucToChar(rand()%4);
		}
		else{
			if(p){
				seq += s1.seq[in];
			}
			else{
				seq += s2.seq[in];
			}
		}

		if(rand()% 15 == 0){
			p = !p;
		}
		++in;
	}

	return sek(seq.c_str());
}

static unsigned char decode(char ch){
	switch(ch){
	case 'A':
       		return 1;
	case 'C':
        	return 2;
    	case 'G':
        	return 4;
    	case 'U':
        	return 8;
	case 'R':
        	return 5;
    	case 'Y':
        	return 10;
    	case 'S':
        	return 6;
    	case 'W':
        	return 9;
    	case 'K':
        	return 12;
    	case 'M':
        	return 3;
    	case 'B':
        	return 14;
    	case 'D':
        	return 13;
    	case 'H':
        	return 11;
    	case 'V':
        	return 7;
    	case 'N':
        	return 15;
	case '(':
		return 16 + 15;
	case '<':
		return 16;
	case ')':
		return 32 + 15;
	case '>':
		return 32;
    	case '.':
		return 64 + 15;
	case '@':
		return 127;
	case '%':
		return 64 + 32 + 16;
	case '#':
		return 32 + 16;
	default:
		return 0;
	}
}

static bool loadTemplate(const char* t){
	if(t[0] == '^'){
		++t;
		begi = true;
	}
	int i = -1;

	while(t[++i]){
		unsigned char res = decode(toupper(t[i]));
		if(!res){
			if(t[i] == '$' && t[i+1] == 0){
				endi = true;
				break;
			}
			return false;
		}

		if(res < 16){ //nucleotide
			temp.push_back(res + 64);
		}
		else{
			if(t[i] == '<' || t[i] == '>' || t[i] == '%' || t[i] == '#'){
				unsigned char res2 = decode(toupper(t[i + 1]));
				if(!res2 || res2 > 15){
					return false;
				}
				temp.push_back(res | res2);
				++i;
			}
			else{
				temp.push_back(res);
			}
		}


	}

	return true;
}

int main(int argc, char** argv){
	puts("Program do wyszukiwania sekwencji RNA z wykorzystaniem algorytmu genetycznego");
	puts("Autor: Jaroslaw Synak");
	//puts("2020-01-29");
	puts("Program wykorzystuje ViennaRNA Package - https://www.tbi.univie.ac.at/RNA/\n\n");




	if(argc < 3){
		fprintf(stderr, "Too few arguments!\n");
		fprintf(stderr, "./RNASearchGen <template> <length>\n");
		return 1;
	}

	if(!loadTemplate(argv[1])){
		fprintf(stderr, "Bad template!\n");
		return 1;
	}


	unsigned int s = atoi(argv[2]);

	if(s > 2048){
		fprintf(stderr, "Too big sequence length! Current limit is 2048!\n");
		return 1;
	}

	#pragma omp parallel
	{

		srand(time(0) ^ omp_get_thread_num()); rand();
		#pragma omp for	
		for(int i = 0; i < M; ++i){
			curr[i] = sek(s);
		}
	}

	sort(curr, curr + M);

	sek best = curr[0];
	int besti = 0;


	printf("%s %i %f gen=0\n", curr[0].seq.c_str(), curr[0].result, curr[0].mfe);

	for(int i = 1; true; ++i){
		if(i - besti >= 300){
			break;
		}
		if(i - besti == 100){
			mutation = 100;
			b = 1500;
		}
		else if(i - besti == 200){
			mutation = 256;
			b = 100;
		}
		#pragma omp parallel
		{
			srand(time(0) ^ omp_get_thread_num()); rand();
			#pragma omp for
			for(int j = 0; j < M; ++j){
				int in1 = choose();
				int in2 = choose();

				curr2[j] = cross(curr[in1], curr[in2]);
			}
		}
		sort(curr2, curr2 + M);
		if(curr2[0] < best){
			mutation = 256;
			b = 100;
			best = curr2[0];
			besti = i;
			printf("%s %i %f gen=%i\n", curr2[0].seq.c_str(), curr2[0].result, curr2[0].mfe, i);
		}

		swap(curr, curr2);
	}



	return 0;
}
