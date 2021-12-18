#include <bits/stdc++.h>
#include <random>
#ifdef _MSC_VER
#include <ppl.h>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>
#else
#pragma GCC target("avx2")
#pragma GCC optimize("O3")
#pragma GCC optimize("unroll-loops")
#endif

/** compro_io **/

/* tuple */
// out
namespace aux {
    template<typename T, unsigned N, unsigned L>
    struct tp {
        static void output(std::ostream& os, const T& v) {
            os << std::get<N>(v) << ", ";
            tp<T, N + 1, L>::output(os, v);
        }
    };
    template<typename T, unsigned N>
    struct tp<T, N, N> {
        static void output(std::ostream& os, const T& v) { os << std::get<N>(v); }
    };
}
template<typename... Ts>
std::ostream& operator<<(std::ostream& os, const std::tuple<Ts...>& t) {
    os << '[';
    aux::tp<std::tuple<Ts...>, 0, sizeof...(Ts) - 1>::output(os, t);
    return os << ']';
}

template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x);

/* pair */
// out
template<class S, class T>
std::ostream& operator<<(std::ostream& os, const std::pair<S, T>& p) {
    return os << "[" << p.first << ", " << p.second << "]";
}
// in
template<class S, class T>
std::istream& operator>>(std::istream& is, const std::pair<S, T>& p) {
    return is >> p.first >> p.second;
}

/* container */
// out
template<class Ch, class Tr, class Container>
std::basic_ostream<Ch, Tr>& operator<<(std::basic_ostream<Ch, Tr>& os, const Container& x) {
    bool f = true;
    os << "[";
    for (auto& y : x) {
        os << (f ? "" : ", ") << y;
        f = false;
    }
    return os << "]";
}
// in
template <
    class T,
    class = decltype(std::begin(std::declval<T&>())),
    class = typename std::enable_if<!std::is_same<T, std::string>::value>::type
>
std::istream& operator>>(std::istream& is, T& a) {
    for (auto& x : a) is >> x;
    return is;
}

/* struct */
template<typename T>
auto operator<<(std::ostream& out, const T& t) -> decltype(out << t.stringify()) {
    out << t.stringify();
    return out;
}

/* setup */
struct IOSetup {
    IOSetup(bool f) {
        if (f) { std::cin.tie(nullptr); std::ios::sync_with_stdio(false); }
        std::cout << std::fixed << std::setprecision(15);
    }
} iosetup(true);

/** string formatter **/
template<typename... Ts>
std::string format(const std::string& f, Ts... t) {
    size_t l = std::snprintf(nullptr, 0, f.c_str(), t...);
    std::vector<char> b(l + 1);
    std::snprintf(&b[0], l + 1, f.c_str(), t...);
    return std::string(&b[0], &b[0] + l);
}

template<typename T>
std::string stringify(const T& x) {
    std::ostringstream oss;
    oss << x;
    return oss.str();
}

/* dump */
#define ENABLE_DUMP
#ifdef ENABLE_DUMP
#define DUMPOUT std::cerr
std::ostringstream DUMPBUF;
#define dump(...) do{DUMPBUF<<"  ";DUMPBUF<<#__VA_ARGS__<<" :[DUMP - "<<__LINE__<<":"<<__FUNCTION__<<"]"<<std::endl;DUMPBUF<<"    ";dump_func(__VA_ARGS__);DUMPOUT<<DUMPBUF.str();DUMPBUF.str("");DUMPBUF.clear();}while(0);
void dump_func() { DUMPBUF << std::endl; }
template <class Head, class... Tail> void dump_func(Head&& head, Tail&&... tail) { DUMPBUF << head; if (sizeof...(Tail) == 0) { DUMPBUF << " "; } else { DUMPBUF << ", "; } dump_func(std::move(tail)...); }
#else
#define dump(...) void(0);
#endif

/* timer */
class Timer {
    double t = 0, paused = 0, tmp;
public:
    Timer() { reset(); }
    static double time() {
#ifdef _MSC_VER
        return __rdtsc() / 3.0e9;
#else
        unsigned long long a, d;
        __asm__ volatile("rdtsc"
            : "=a"(a), "=d"(d));
        return (d << 32 | a) / 3.0e9;
#endif
    }
    void reset() { t = time(); }
    void pause() { tmp = time(); }
    void restart() { paused += time() - tmp; }
    double elapsed_ms() { return (time() - t - paused) * 1000.0; }
} timer;

/* rand */
struct Xorshift {
    uint64_t x = 88172645463325252LL;
    void set_seed(unsigned seed, int rep = 100) { x = uint64_t((seed + 1) * 10007); for (int i = 0; i < rep; i++) next_int(); }
    unsigned next_int() { x = x ^ (x << 7); return x = x ^ (x >> 9); }
    unsigned next_int(unsigned mod) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % mod; }
    unsigned next_int(unsigned l, unsigned r) { x = x ^ (x << 7); x = x ^ (x >> 9); return x % (r - l + 1) + l; } // inclusive
    double next_double() { return double(next_int()) / UINT_MAX; }
} rnd;

/* shuffle */
template<typename T>
void shuffle_vector(std::vector<T>& v, Xorshift& rnd) {
    int n = v.size();
    for (int i = n - 1; i >= 1; i--) {
        int r = rnd.next_int(i);
        std::swap(v[i], v[r]);
    }
}

/* split */
std::vector<std::string> split(std::string str, const std::string& delim) {
    for (char& c : str) if (delim.find(c) != std::string::npos) c = ' ';
    std::istringstream iss(str);
    std::vector<std::string> parsed;
    std::string buf;
    while (iss >> buf) parsed.push_back(buf);
    return parsed;
}

template<typename A, size_t N, typename T> inline void Fill(A(&array)[N], const T& val) {
    std::fill((T*)array, (T*)(array + N), val);
}

template<typename T> bool chmax(T& a, const T& b) { if (a < b) { a = b; return true; } return false; }
template<typename T> bool chmin(T& a, const T& b) { if (a > b) { a = b; return true; } return false; }



using namespace std;

constexpr int N = 20;
constexpr int NN = N * N;

constexpr int di[] = { 0, -1, 0, 1 };
constexpr int dj[] = { 1, 0, -1, 0 };

string rot_char = "lLrR";

//vector<string> rot_str({
//    "l", "L", "r", "R"
//    });

struct TestCase;
using TestCasePtr = shared_ptr<TestCase>;
struct TestCase {
    int si, sj;
    vector<string> h;
    vector<string> v;
    TestCase(istream& in) {
        in >> si >> sj;
        h.resize(N, string(N + 1, '1'));
        v.resize(N + 1, string(N, '1'));
        for (int y = 0; y < N; y++) {
            in >> h[y];
            h[y] = '1' + h[y] + '1';
        }
        for (int y = 1; y < N; y++) {
            in >> v[y];
        }
    }
    string stringify() const {
        ostringstream oss;
        oss << si << ' ' << sj;
        for (const auto& s : h) oss << '\n' << s;
        for (const auto& s : v) oss << '\n' << s;
        return oss.str();
    }
};

string units2str(const string& units) {
    string ret;
    for (char c : units) {
        ret += c;
        ret += 'F';
    }
    return ret;
}

struct Block {
    int rep;
    string units;
    string stringify() const {
        if (rep == 1) return units2str(units);
        else return to_string(rep) + "(" + units2str(units) + ")";
    }
    static Block create_random(Xorshift& rnd) {
        Block block;
        block.rep = 3;
        for (int i = 0; i < 5; i++) {
            block.units.push_back(rot_char[rnd.next_int(4)]);
        }
        return block;
    }
};

string stringify(const vector<Block>& blocks) {
    string res;
    for (const auto& block : blocks) res += block.stringify();
    return res;
}

struct SuperBlock {
    int rep;
    vector<Block> blocks;
    string stringify() const {
        if (rep == 1) return ::stringify(blocks);
        else return to_string(rep) + "(" + ::stringify(blocks) + ")";
    }
    static SuperBlock create_random(Xorshift& rnd) {
        SuperBlock sblock;
        sblock.rep = 9;
        for (int i = 0; i < 2; i++) {
            sblock.blocks.push_back(Block::create_random(rnd));
        }
        return sblock;
    }
};

struct State {
    TestCasePtr tc;
    int turn;
    int dir, i, j;
    int num_seen;
    bool seen[N][N];
    State(TestCasePtr tc) : tc(tc) { reset(); }
    void reset() {
        turn = 0;
        dir = 1;
        i = tc->si;
        j = tc->sj;
        memset(seen, 0, sizeof(bool) * NN);
        seen[i][j] = true;
        num_seen = 1;
    }
    bool is_blocked() const {
        return
            (dir == 0 && tc->h[i][j + 1] == '1') ||
            (dir == 1 && tc->v[i][j] == '1') ||
            (dir == 2 && tc->h[i][j] == '1') ||
            (dir == 3 && tc->v[i + 1][j] == '1');
    }
    void proceed(const char c) {
        bool blocked = is_blocked();
        if (c == 'L') {
            dir = (dir + 1) & 3;
        }
        else if (c == 'R') {
            dir = (dir + 3) & 3;
        }
        else if (c == 'l') {
            if (blocked) dir = (dir + 1) & 3;
        }
        else if (c == 'r') {
            if (blocked) dir = (dir + 3) & 3;
        }
        else if (c == 'F') {
            if (!blocked) {
                i += di[dir];
                j += dj[dir];
                if (!seen[i][j]) {
                    seen[i][j] = true;
                    num_seen++;
                }
            }
        }
        turn++;
    }
    void proceed(const Block& b) {
        for (int i = 0; i < b.rep; i++) {
            for (char unit : b.units) {
                proceed(unit);
                proceed('F');
            }
        }
    }
    void proceed(const SuperBlock& sb) {
        for (int i = 0; i < sb.rep; i++) {
            for (const auto& b : sb.blocks) {
                proceed(b);
            }
        }
    }
};

string decode(const string& src) {
    string dst;
    for (int i = 0; i < src.size(); i++) {
        if (isdigit(src[i])) {
            dst += string(src[i] - '0', src[i + 1]);
            i++;
        }
        else {
            dst += src[i];
        }
    }
    return dst;
}

int solve(istream& in, ostream& out, bool no_output = false) {
    Xorshift rnd;
    Timer timer;
    using pii = pair<int, int>;

    auto tc = make_shared<TestCase>(in);

    State best_state(tc);
    SuperBlock best_sblock;
    State prev_state(best_state);
    string ans;

    while (best_state.num_seen < NN && best_state.turn < 5000 && timer.elapsed_ms() < 1800) {
        for (int t = 0; t < 50000; t++) {
            auto sblock = SuperBlock::create_random(rnd);
            State state(prev_state);
            state.proceed(sblock);
            if (chmax(best_state.num_seen, state.num_seen)) {
                best_state = state;
                best_sblock = sblock;
            }
        }
        prev_state = best_state;
        ans += best_sblock.stringify();
    }

    if (!no_output) out << ans << endl;

    int score = ((best_state.num_seen == NN && best_state.turn <= 5000) ? best_state.num_seen + (int)round(1e8 / (100.0 + ans.size())) : best_state.num_seen);

    dump(timer.elapsed_ms(), score);
    return score;
}

#ifdef _MSC_VER
long long batch_test(int num_seeds) {

    vector<long long> scores(num_seeds, 0);

    concurrency::parallel_for(0, num_seeds, [&](int seed) {
        string in_file = format("C:\\dev\\heuristic\\tasks\\HTTF2022final\\tools\\in\\%04d.txt", seed);
        string out_file = format("C:\\dev\\heuristic\\tasks\\HTTF2022final\\tools\\out\\%04d.txt", seed);
        ifstream ifs(in_file);
        ofstream ofs(out_file);
        int score = solve(ifs, ofs, false);
        scores[seed] = score;
        }, concurrency::simple_partitioner(num_seeds / 5));

    long long score_sum = accumulate(scores.begin(), scores.end(), 0LL);
    dump(score_sum);
    return score_sum;
}

int single_test() {
    string in_file = "C:\\dev\\heuristic\\tasks\\HTTF2022final\\tools\\in\\0000.txt";
    string out_file = "C:\\dev\\heuristic\\tasks\\HTTF2022final\\tools\\out\\0000.txt";
    ifstream ifs(in_file);
    ofstream ofs(out_file);
    int score = solve(ifs, ofs);
    dump(score);
    return score;
}

int main() {
    //batch_test(100);
    single_test();
}

#else

int main() {

    std::istream& in = std::cin;
    std::ostream& out = std::cout;

    solve(in, out);

    return 0;
}

#endif