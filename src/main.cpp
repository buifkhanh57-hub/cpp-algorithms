// cpp-algorithms — classic algorithms and a BST in modern C++.
// Build: make && ./build/algorithms

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <memory>
#include <optional>
#include <queue>
#include <string>
#include <vector>

namespace algo {

// ------------------------------------------------------------ sorting ------

template <typename T>
void quick_sort(std::vector<T>& a, int lo, int hi) {
    if (lo >= hi) return;
    const T pivot = a[lo + (hi - lo) / 2];
    int i = lo, j = hi;
    while (i <= j) {
        while (a[i] < pivot) i++;
        while (a[j] > pivot) j--;
        if (i <= j) std::swap(a[i++], a[j--]);
    }
    quick_sort(a, lo, j);
    quick_sort(a, i, hi);
}

template <typename T>
void quick_sort(std::vector<T>& a) {
    quick_sort(a, 0, static_cast<int>(a.size()) - 1);
}

std::vector<int> merge_sorted(const std::vector<int>& a, const std::vector<int>& b) {
    std::vector<int> out;
    out.reserve(a.size() + b.size());
    std::merge(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(out));
    return out;
}

// ----------------------------------------------------------- searching -----

template <typename T>
int binary_search_idx(const std::vector<T>& a, const T& target) {
    int lo = 0, hi = static_cast<int>(a.size()) - 1;
    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] == target) return mid;
        if (a[mid] < target) lo = mid + 1; else hi = mid - 1;
    }
    return -1;
}

template <typename T>
int lower_bound_idx(const std::vector<T>& a, const T& target) {
    int lo = 0, hi = static_cast<int>(a.size());
    while (lo < hi) {
        int mid = lo + (hi - lo) / 2;
        if (a[mid] < target) lo = mid + 1; else hi = mid;
    }
    return lo;
}

// ---------------------------------------------------------------- BST ------

class Bst {
public:
    void insert(int64_t key) { root_ = insertRec(std::move(root_), key); }

    bool contains(int64_t key) const {
        const Node* cur = root_.get();
        while (cur) {
            if (key == cur->key) return true;
            cur = key < cur->key ? cur->left.get() : cur->right.get();
        }
        return false;
    }

    std::vector<int64_t> in_order() const {
        std::vector<int64_t> out;
        inOrderRec(root_.get(), out);
        return out;
    }

    std::optional<int64_t> min() const {
        if (!root_) return std::nullopt;
        const Node* cur = root_.get();
        while (cur->left) cur = cur->left.get();
        return cur->key;
    }

private:
    struct Node {
        int64_t key;
        std::unique_ptr<Node> left, right;
        explicit Node(int64_t k) : key(k) {}
    };

    std::unique_ptr<Node> insertRec(std::unique_ptr<Node> node, int64_t key) {
        if (!node) return std::make_unique<Node>(key);
        if (key < node->key) node->left = insertRec(std::move(node->left), key);
        else if (key > node->key) node->right = insertRec(std::move(node->right), key);
        return node; // ignore duplicates
    }

    static void inOrderRec(const Node* node, std::vector<int64_t>& out) {
        if (!node) return;
        inOrderRec(node->left.get(), out);
        out.push_back(node->key);
        inOrderRec(node->right.get(), out);
    }

    std::unique_ptr<Node> root_;
};

// --------------------------------------------------------- misc problems ---

long long max_subarray_sum(const std::vector<int>& a) {
    // Kadane's algorithm
    long long best = 0, cur = 0;
    for (int x : a) {
        cur = std::max<long long>(cur + x, x);
        best = std::max(best, cur);
    }
    return best;
}

std::string reverse_words(const std::string& s) {
    std::vector<std::string> words;
    std::string cur;
    for (char c : s) {
        if (std::isspace(static_cast<unsigned char>(c))) {
            if (!cur.empty()) words.push_back(cur);
            cur.clear();
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) words.push_back(cur);
    std::string out;
    for (auto it = words.rbegin(); it != words.rend(); ++it) {
        if (it != words.rbegin()) out += ' ';
        out += *it;
    }
    return out;
}

int count_islands(std::vector<std::vector<int>> grid) {
    // Flood fill over a copy; counts connected components of 1s.
    if (grid.empty()) return 0;
    const int rows = static_cast<int>(grid.size());
    const int cols = static_cast<int>(grid[0].size());
    int islands = 0;
    for (int r = 0; r < rows; ++r) {
        for (int c = 0; c < cols; ++c) {
            if (grid[r][c] != 1) continue;
            ++islands;
            std::queue<std::pair<int, int>> q;
            q.push({r, c});
            grid[r][c] = 0;
            while (!q.empty()) {
                auto [cr, cc] = q.front(); q.pop();
                const int dr[] = {-1, 1, 0, 0};
                const int dc[] = {0, 0, -1, 1};
                for (int d = 0; d < 4; ++d) {
                    int nr = cr + dr[d], nc = cc + dc[d];
                    if (nr < 0 || nr >= rows || nc < 0 || nc >= cols) continue;
                    if (grid[nr][nc] == 1) {
                        grid[nr][nc] = 0;
                        q.push({nr, nc});
                    }
                }
            }
        }
    }
    return islands;
}

} // namespace algo

int main() {
    std::vector<int> data{42, -7, 13, 0, 99, -100, 55, 8, 23, -4};
    algo::quick_sort(data);
    assert(std::is_sorted(data.begin(), data.end()));
    assert(algo::binary_search_idx(data, 23) != -1);
    assert(algo::binary_search_idx(data, 24) == -1);
    assert(algo::lower_bound_idx(data, 0) == 3);

    auto merged = algo::merge_sorted({1, 3, 5}, {2, 4, 6});
    assert((merged == std::vector<int>{1, 2, 3, 4, 5, 6}));

    algo::Bst tree;
    for (int64_t k : {50, 30, 70, 20, 40, 60, 80}) tree.insert(k);
    assert(tree.contains(40) && !tree.contains(35));
    assert(tree.min().value_or(-1) == 20);
    auto sorted_out = tree.in_order();
    assert(std::is_sorted(sorted_out.begin(), sorted_out.end()));

    assert(algo::max_subarray_sum({-2, 1, -3, 4, -1, 2, 1, -5, 4}) == 6);
    assert(algo::reverse_words("  the quick   brown fox ") == "fox brown quick the");
    assert(algo::count_islands({{1, 1, 0}, {0, 1, 0}, {0, 0, 1}}) == 2);

    std::cout << "All C++ algorithm checks passed.\n";
    for (int64_t v : sorted_out) std::cout << v << ' ';
    std::cout << '\n';
    return 0;
}
