#include "inversionmesh.h"

void InversionMesh::clear() {
    elements_.clear();
}

void InversionMesh::addElement(const Element& e) {
    elements_.push_back(e);
}

size_t InversionMesh::size() const {
    return elements_.size();
}

const InversionMesh::Element& InversionMesh::element(size_t i) const {
    return elements_[i];
}

InversionMesh::Element& InversionMesh::element(size_t i) {
    return elements_[i];
}

const std::vector<InversionMesh::Element>& InversionMesh::elements() const {
    return elements_;
}
