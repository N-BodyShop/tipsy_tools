/*
 ** Header record definitions for Hugh's code.
 */
#define REC_HEADER1		24
#define REC_HEADER2		400

struct header1 {
	int iRecStart;
	int nobj;
	int ips;
	float aexp;
	float p;
	float dpnp1;
	float dpn;
	int iRecEnd;
	};

struct header2 {
	int iRecStart;
	union ibuf_union {
		struct common_block {
			int nobj;
			int ips;
			int istart;
			int irun;
			int ipout;
			int ipstop;
			int ipdump;
			int iseed1;
			int iseed2;
			int L;
			int iru;
			int nlmx;
			float alpha;
			float perr;
			float aexp;
			float p;
			float dpnp1;
			float dpn;
			float Ust;
			float Tst;
			float Uh;
			float Usum;
			float Usm;
			float T;
			float Tsum;
			float delI;
			float delC;
			float dfCU;
			float rLbox;
			float delta;
			float soft;
			float sftinit;
			float sftmin;
			} common;
		int ibuf[100];
		} u;
	int iRecEnd;
	};

