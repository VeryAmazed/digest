#include <algorithm>
#include <cctype>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <deque>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <list>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;
typedef long long ll;
// if you end up using long double, you need to set the floating point notation
// to fixed, and set the percision to be very high
typedef long double ld;

// contrsuct umaps like this, unordered_map<long long, int, custom_hash>
// safe_map; FIXED_RANDOM is static so it doesn not get redeclared between
// function calls
struct custom_hash {
	static uint64_t splitmix64(uint64_t x) {
		// http://xorshift.di.unimi.it/splitmix64.c
		x += 0x9e3779b97f4a7c15;
		x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9;
		x = (x ^ (x >> 27)) * 0x94d049bb133111eb;
		return x ^ (x >> 31);
	}

	size_t operator()(uint64_t x) const {

		static const uint64_t FIXED_RANDOM =
			chrono::steady_clock::now().time_since_epoch().count();
		return splitmix64(x + FIXED_RANDOM);
	}
};

#define INF 2001001001
#define INF2 2e18
#define MOD 1000000007

#define f0r(a, b) for (long long a = 0; a < b; a++)
#define f1r(a, b, c) for (long long a = b; a < c; a++)
#define max3(a, b, c) max(a, max(b, c))
#define min3(a, b, c) min(a, min(b, c))
#define pb push_back
#define pf push_front
#define f first
#define s second
#define mp make_pair
#define pll pair<ll, ll>
#define pii pair<int, int>
#define tp make_tuple

// first four are north, west, east ,south
int dir1[] = {1, 0, -1, 0, 1, 1, -1, -1};
int dir2[] = {0, 1, 0, -1, 1, -1, 1, -1};

int main() {

	// use this if you read in from a file

	// freopen("in.txt", "r", stdin);
	freopen("ACTG.txt", "w", stdout);

	char chars[4] = {'A', 'C', 'T', 'G'};
	for (int i = 0; i < 100; i++) {
		string temp;
		for (int j = 0; j < 1e5; j++) {
			int curr = rand() % 4;
			temp.pb(chars[curr]);
		}
		cout << temp << endl;
		cout << endl;
	}
	freopen("non-ACTG.txt", "w", stdout);
	char chars2[5] = {'A', 'C', 'T', 'G', 'N'};
	for (int i = 0; i < 100; i++) {
		string temp;
		for (int j = 0; j < 1e5; j++) {
			int curr = rand() % 33;
			if (curr == 32) {
				temp.pb(chars2[4]);
			} else {
				curr %= 4;
				temp.pb(chars2[curr]);
			}
		}
		cout << temp << endl;
		cout << endl;
	}
	return 0;
}