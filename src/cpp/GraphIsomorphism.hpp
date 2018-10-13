#ifndef GRAPH_ISOMORPHISM_HPP
#define GRAPH_ISOMORPHISM_HPP

#include <nauty.h>
#include <memory>
#include <vector>

template <typename T>
using Vector = std::vector<T>;

class GraphIsomorphism {
public:
	GraphIsomorphism() : _isoLabel(0), _qty(1) {}
	GraphIsomorphism(std::unique_ptr<graph[]> canonicalLabel) : _isoLabel(0), _qty(1), _canonicalLabel(std::move(canonicalLabel)) {}

	inline void set_iso_label(int isoLabel) {
		_isoLabel = isoLabel;
	}

	inline int get_iso_label() const {
		return _isoLabel;
	}

	inline void add_qty(int qty) {
		_qty += qty;
	}

	inline int get_qty() const {
		return _qty;
	}

	inline graph * get_canonical_label() const {
		return _canonicalLabel.get();
	}
private:
	int _isoLabel;
	int _qty;
	std::unique_ptr<graph[]> _canonicalLabel;
};

#endif /* GRAPH_ISOMORPHISM_HPP */
