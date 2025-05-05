#include <benchmark/benchmark.h>
#include <digest/data_structure.hpp>

#include <array>
#include <cstdint>
#include <iostream>
#include <random>

// segtee wins at 12
// naive2 wins at 17

namespace digest::ds {

// sshash implementation of monoqueue (ours is better)
//template <typename T>
//struct fixed_size_deque {
//	fixed_size_deque(uint32_t size) : m_begin(0), m_end(0), m_size(0), m_buffer(size) {
//	}
//
//	bool empty() const { return m_size == 0; }
//
//	void push_back(T const& val) {
//		m_buffer[m_end] = val;
//		if (++m_end == m_buffer.size()) m_end = 0;
//		++m_size;
//	}
//
//	void pop_front() {
//		if (++m_begin == m_buffer.size()) m_begin = 0;
//		--m_size;
//	}
//
//	void pop_back() {
//		if (m_end == 0) m_end = m_buffer.size();
//		--m_end;
//		--m_size;
//	}
//
//	T const& front() const { return m_buffer[m_begin]; }
//
//	T const& back() const {
//		if (m_end == 0) return m_buffer.back();
//		return m_buffer[m_end - 1];
//	}
//
//private:
//	uint32_t m_begin;
//	uint32_t m_end;
//	uint32_t m_size;
//	std::vector<T> m_buffer;
//};
//
//struct minimizer_enumerator {
//
//	minimizer_enumerator(uint32_t k)
//		: m_k(k)
//		, m_position(0)
//		, m_q(k)
//	{}
//
//	uint32_t min() { return m_q.front().value; }
//
//	void insert(uint32_t pos, uint32_t hash) {
//		/* Removes from front elements which are no longer in the window */
//		while (!m_q.empty() and m_position >= m_k and
//			   m_q.front().position <= m_position  - m_k) {
//			m_q.pop_front();
//		}
//		/* Removes from back elements which are no longer useful */
//		while (!m_q.empty() and hash < m_q.back().hash) m_q.pop_back();
//
//		m_q.push_back({hash, m_position, pos});
//		m_position += 1;
//	}
//
//private:
//	uint32_t m_k;
//	uint32_t m_position;
//
//	struct mmer_t {
//		mmer_t() {}
//		mmer_t(uint32_t hash, uint32_t position, uint32_t value)
//			: hash(hash), position(position), value(value) {}
//		uint32_t hash, position, value;
//	};
//
//	/* NOTE: we could use a std::deque<mmer_t> here,
//			 but std::deque is terribly space-inefficient. */
//	fixed_size_deque<mmer_t> m_q;
//};

template <uint32_t k> struct MonoQueue {
	int head = 0, tail = 0;
	// one extra slot so empty() works when full
	// {hash, index, time}
	std::array<std::array<uint32_t, k + 1>, 3> queue;
	uint32_t time = 0;

	bool empty() { return head == tail; }

	MonoQueue(uint32_t) {
		queue[2].fill(0); // necessary to avoid 1 in a billion collision
		queue[0].fill(0); // gets rid of warning
	}
	MonoQueue(const MonoQueue &other) = default;
	MonoQueue &operator=(const MonoQueue &other) = default;

	void insert(uint32_t index, uint32_t hash) {
		if (queue[2][head] == time - k) {
			if (++head == k + 1) {
				head = 0;
			}
		}

		while (not empty() and queue[0][tail == 0 ? k : tail - 1] >= hash) {
			if (--tail == -1) {
				tail = k;
			}
		}

		queue[0][tail] = hash;
		queue[1][tail] = index;
		queue[2][tail] = time++;

		if (++tail == k + 1)
			tail = 0;
	}

	uint32_t min() { return queue[1][head]; }

	// void min_syncmer(std::vector<uint32_t> &vec) {
	// 	if (queue[2][head] == time - k or queue[2][head] == time - 1) {
	// 		vec.emplace_back(queue[1][head]);
	// 	}
	// }
};

template <uint32_t k> struct Set {
	std::set<uint64_t> mset;
	std::array<std::set<uint64_t>::iterator, k> vec;
	int i = 0;

	Set(uint32_t) {
		// edge case where hash == all 1's can cause a set collision.
		for (unsigned i = 0; i < k; i++) {
			vec[i] = mset.emplace(i).first;
		}
	}
	Set(const Set &other) = delete; // have to copy over iterators
	Set &operator=(const Set &other) = delete;

	void insert(uint32_t index, uint32_t hash) {
		mset.erase(vec[i]);

		vec[i] = mset.emplace((uint64_t)~hash << 32 | index).first;
		if (++i == k)
			i = 0;
	}

	uint32_t min() { return *mset.rbegin(); }
};

} // namespace digest::ds

const int INPUT_SIZE = 1e7;
std::array<uint32_t, 2 * INPUT_SIZE> hashes;

std::map<int, std::map<int, std::array<uint32_t, INPUT_SIZE>>> all;

void setupInput() {
	std::random_device rd;	// seed
	std::mt19937 gen(rd()); // generator
	std::uniform_int_distribution<uint32_t> distrib(0,
													UINT32_MAX); // [0, 2**32]
	for (uint32_t &h : hashes) {
		h = distrib(gen);
	}

	// edge test for ties
	// for (int i = 0; i < 2*INPUT_SIZE; i++) {
	// 	hashes[i] = hashes[0];
	// }
}

template <int k, class T, int out> static void BM(benchmark::State &state) {
	auto &temp = all[out][k];
	for (auto _ : state) {
		T ds(k);
		for (int i = 0; i < k - 1; i++) {
			ds.insert(i, hashes[i]);
		}
		for (int i = 0; i < INPUT_SIZE; i++) {
			ds.insert(i + k - 1, hashes[i + k - 1]);
			temp[i] = ds.min();
		}
		benchmark::ClobberMemory();
	}
}

#define test(name, out)                                                        \
	BENCHMARK_TEMPLATE(BM, 4, name<4>, out);                                   \
	BENCHMARK_TEMPLATE(BM, 5, name<5>, out);                                   \
	BENCHMARK_TEMPLATE(BM, 8, name<8>, out);                                   \
	BENCHMARK_TEMPLATE(BM, 9, name<9>, out);                                   \
	BENCHMARK_TEMPLATE(BM, 12, name<12>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 16, name<16>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 17, name<17>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 32, name<32>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 33, name<33>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 64, name<64>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 96, name<96>, out);                                 \
	BENCHMARK_TEMPLATE(BM, 128, name<128>, out);                               \
	BENCHMARK_TEMPLATE(BM, 256, name<256>, out);                               \
	BENCHMARK_TEMPLATE(BM, 512, name<512>, out);                               \
	BENCHMARK_TEMPLATE(BM, 1024, name<1024>, out);

#define test2(name, out)                                                       \
	BENCHMARK_TEMPLATE(BM, 4, name, out);                                      \
	BENCHMARK_TEMPLATE(BM, 5, name, out);                                      \
	BENCHMARK_TEMPLATE(BM, 8, name, out);                                      \
	BENCHMARK_TEMPLATE(BM, 9, name, out);                                      \
	BENCHMARK_TEMPLATE(BM, 12, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 16, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 17, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 32, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 33, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 64, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 96, name, out);                                     \
	BENCHMARK_TEMPLATE(BM, 128, name, out);                                    \
	BENCHMARK_TEMPLATE(BM, 256, name, out);                                    \
	BENCHMARK_TEMPLATE(BM, 512, name, out);                                    \
	BENCHMARK_TEMPLATE(BM, 1024, name, out);

test(digest::ds::Naive, 0) test(digest::ds::Naive2, 1)
	test(digest::ds::MonoQueue, 2) test(digest::ds::SegmentTree, 3)
		test(digest::ds::Set, 4) test2(digest::ds::Adaptive, 5)
			test2(digest::ds::Adaptive64, 6)
				int main(int argc, char **argv) {
	setupInput();
	benchmark::Initialize(&argc, argv);
	benchmark::RunSpecifiedBenchmarks();

	// sanity check
	for (auto &[_, m] : all) {
		assert(m == all.begin()->second);
	}

	std::cout << "Passed Asserts!" << std::endl;
}
