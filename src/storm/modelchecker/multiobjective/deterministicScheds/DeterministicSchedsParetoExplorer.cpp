#include <sstream>
#include <algorithm>


#include "storm/modelchecker/multiobjective/deterministicScheds/DeterministicSchedsParetoExplorer.h"
#include "storm/storage/geometry/coordinates.h"
#include "storm/models/sparse/MarkovAutomaton.h"
#include "storm/models/sparse/Mdp.h"
#include "storm/models/sparse/StandardRewardModel.h"
#include "storm/modelchecker/multiobjective/MultiObjectivePostprocessing.h"
#include "storm/modelchecker/results/ExplicitParetoCurveCheckResult.h"
#include "storm/environment/modelchecker/MultiObjectiveModelCheckerEnvironment.h"

#include "storm/utility/export.h"
#include "storm/utility/solver.h"

#include "storm/exceptions/UnexpectedException.h"
#include "storm/exceptions/InvalidOperationException.h"


namespace storm {
    namespace modelchecker {
        namespace multiobjective {
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::Point(std::vector<GeometryValueType> const& coordinates) : coordinates(coordinates), paretoOptimal(false), onFacet(false) {
                STORM_LOG_ASSERT(!this->coordinates.empty(), "Points with dimension 0 are not supported");
            }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::Point(std::vector<GeometryValueType>&& coordinates) : coordinates(std::move(coordinates)), paretoOptimal(false), onFacet(false) {
                STORM_LOG_ASSERT(!this->coordinates.empty(), "Points with dimension 0 are not supported");
            }
            
            template <class SparseModelType, typename GeometryValueType>
            std::vector<GeometryValueType>& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::get() {
                return coordinates;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            std::vector<GeometryValueType> const& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::get() const {
                return coordinates;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            uint64_t DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::dimension() const {
                STORM_LOG_ASSERT(!coordinates.empty(), "Points with dimension 0 are not supported");
                return coordinates.size();
            }
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::DominanceResult DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::getDominance(Point const& other) const {
                STORM_LOG_ASSERT(this->dimension() == other.dimension(), "Non-Equal dimensions of points: [" << this->toString() << "] vs. [" << other.toString() << "]");
                auto thisIt = this->get().begin();
                auto otherIt = other.get().begin();
                auto thisItE = this->get().end();
                
                // Find the first entry where the points differ
                while (*thisIt == *otherIt) {
                    ++thisIt;
                    ++otherIt;
                    if (thisIt == thisItE) {
                        return DominanceResult::Equal;
                    }
                }
                    
                if (*thisIt > *otherIt) {
                    // *this might dominate other
                    for (++thisIt, ++otherIt; thisIt != thisItE; ++thisIt, ++otherIt) {
                        if (*thisIt < *otherIt) {
                            return DominanceResult::Incomparable;
                        }
                    }
                    return DominanceResult::Dominates;
                } else {
                    assert(*thisIt < *otherIt);
                    // *this might be dominated by other
                    for (++thisIt, ++otherIt; thisIt != thisItE; ++thisIt, ++otherIt) {
                        if (*thisIt > *otherIt) {
                            return DominanceResult::Incomparable;
                        }
                    }
                    return DominanceResult::Dominated;
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::setParetoOptimal(bool value) {
                paretoOptimal = value;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            bool DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::isParetoOptimal() const {
                return paretoOptimal;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::setOnFacet(bool value) {
                onFacet = value;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            bool DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::liesOnFacet() const {
                return onFacet;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            std::string DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point::toString(bool convertToDouble) const {
                std::stringstream out;
                bool first = true;
                for (auto const& pi : this->get()) {
                    if (first) {
                        first = false;
                    } else {
                        out << ", ";
                    }
                    if (convertToDouble) {
                        out << storm::utility::convertNumber<double>(pi);
                    } else {
                        out << pi;
                    }
                }
                return out.str();
            }
            
      //      template <class SparseModelType, typename GeometryValueType>
      //      bool operator<(typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point const& lhs, typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point const& rhs) {
      //          STORM_LOG_ASSERT(lhs.dimension() == rhs.dimension(), "Non-Equal dimensions of points: " << lhs << " vs. " << rhs);
      //          for (uint64_t i = 0; i < lhs.dimension(); ++i) {
      //              if (lhs.get()[i] < rhs.get()[i]) {
      //                  return true;
      //              } else if (lhs.get()[i] != rhs.get()[i]) {
      //                  return false;
      //              }
      //          }
      //          return false;
      //      }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::Pointset() : currId(1) {
                // Intentionally left empty
            }
        
            template <class SparseModelType, typename GeometryValueType>
            boost::optional<typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::PointId> DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::addPoint(Environment const& env, Point&& point) {
            
                // Find dominated and dominating points
                auto pointsIt = points.begin();
                while (pointsIt != points.end()) {
                    switch (point.getDominance(pointsIt->second)) {
                        case Point::DominanceResult::Incomparable:
                            // Nothing to be done for this point
                            ++pointsIt;
                            break;
                        case Point::DominanceResult::Dominates:
                            // Found a point in the set that is dominated by the new point, so we erase it
                            if (pointsIt->second.isParetoOptimal()) {
                                STORM_LOG_WARN("Potential precision issues: Found a point that dominates another point which was flagged as pareto optimal. Distance of points is " << std::sqrt(storm::utility::convertNumber<double>(storm::storage::geometry::squaredEuclideanDistance(pointsIt->second.get(), point.get()))));
                                point.setParetoOptimal(true);
                            }
                            if (pointsIt->second.liesOnFacet()) {
                                // Do not erase points that lie on a facet
                                ++pointsIt;
                            } else {
                                pointsIt = points.erase(pointsIt);
                            }
                            break;
                        case Point::DominanceResult::Dominated:
                            // The new point is dominated by another point.
                            return boost::none;
                        case Point::DominanceResult::Equal:
                            if (point.isParetoOptimal()) {
                                pointsIt->second.setParetoOptimal();
                            }
                            if (point.liesOnFacet()) {
                                pointsIt->second.setOnFacet();
                            }
                            return pointsIt->first;
                    }
                }
                
                if (env.modelchecker().multi().isPrintResultsSet()) {
                    std::cout << "## achievable point: [" << point.toString(true) << "]" << std::endl;
                }
                
                points.emplace_hint(points.end(), currId, std::move(point));
                return currId++;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Point const& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::getPoint(PointId const& id) const {
                return points.at(id);
            }
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::iterator_type DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::begin() const {
                return points.begin();
            }
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::iterator_type DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::end() const {
                return points.end();
            }
            
            template <class SparseModelType, typename GeometryValueType>
            uint64_t DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::size() const {
                return points.size();
            }
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Polytope DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::downwardClosure() const {
                std::vector<std::vector<GeometryValueType>> pointsAsVector;
                pointsAsVector.reserve(size());
                for (auto const& p : points) {
                    pointsAsVector.push_back(p.second.get());
                }
                return storm::storage::geometry::Polytope<GeometryValueType>::createDownwardClosure(std::move(pointsAsVector));
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::collectPointsInPolytope(std::set<PointId>& collectedPoints, Polytope const& polytope) {
                for (auto const& p : points) {
                    if (polytope->contains(p.second.get())) {
                        collectedPoints.insert(p.first);
                    }
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Pointset::printToStream(std::ostream& out, bool includeIDs, bool convertToDouble) {
                for (auto const& p : this->points) {
                    if (includeIDs) {
                        out << p.first << ": [" << p.second.toString(convertToDouble) << "]" << std::endl;
                    } else {
                        out << p.second.toString(convertToDouble) << std::endl;
                    }
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::Facet(storm::storage::geometry::Halfspace<GeometryValueType> const& halfspace) : halfspace(halfspace) {
                // Intentionally left empty
            }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::Facet(storm::storage::geometry::Halfspace<GeometryValueType>&& halfspace) : halfspace(std::move(halfspace)) {
                // Intentionally left empty
            }
            
            template <class SparseModelType, typename GeometryValueType>
            storm::storage::geometry::Halfspace<GeometryValueType> const& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::getHalfspace() const {
                return halfspace;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::addPoint(PointId const& pointId, Point const& point) {
                inducedSimplex = nullptr;
                GeometryValueType product = storm::utility::vector::dotProduct(getHalfspace().normalVector(), point.get());
                if (product != getHalfspace().offset()) {
                    if (product < getHalfspace().offset()) {
                        STORM_LOG_DEBUG("The point on the facet actually has distance " << storm::utility::convertNumber<double>(getHalfspace().euclideanDistance(point.get())));
                    } else {
                        STORM_LOG_DEBUG("Halfspace of facet is shifted by " << storm::utility::convertNumber<double>(getHalfspace().euclideanDistance(point.get())) << " to capture all points that are supposed to lie on the facet.");
                        halfspace.offset() = product;
                    }
                }
                paretoPointsOnFacet.push_back(pointId);
            }
            
            template <class SparseModelType, typename GeometryValueType>
            std::vector<typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::PointId> const& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::getPoints() const {
                return paretoPointsOnFacet;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            uint64_t DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::getNumberOfPoints() const {
                return paretoPointsOnFacet.size();
            }
            
            
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Polytope const& DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Facet::getInducedSimplex(Pointset const& pointset, std::vector<GeometryValueType> const& referenceCoordinates) {
                if (!inducedSimplex) {
                    std::vector<std::vector<GeometryValueType>> vertices = {referenceCoordinates};
                    for (auto const& pId : paretoPointsOnFacet) {
                        vertices.push_back(pointset.getPoint(pId).get());
                    }
                    // This facet might lie at the 'border', which means that the downward closure has to be taken in some directions
                    storm::storage::BitVector dimensionsForDownwardClosure = storm::utility::vector::filterZero(this->halfspace.normalVector());
                    STORM_LOG_ASSERT(dimensionsForDownwardClosure.getNumberOfSetBits() + vertices.size() >= halfspace.normalVector().size() + 1, "The number of points on the facet is insufficient");
                    if (dimensionsForDownwardClosure.empty()) {
                        inducedSimplex = storm::storage::geometry::Polytope<GeometryValueType>::create(vertices);
                    } else {
                        inducedSimplex = storm::storage::geometry::Polytope<GeometryValueType>::createSelectiveDownwardClosure(vertices, dimensionsForDownwardClosure);
                    }
                }
                return inducedSimplex;
            }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::FacetAnalysisContext::FacetAnalysisContext(Facet& f) : facet(f) {
                // Intentionally left empty
            }
            
            template <class SparseModelType, typename GeometryValueType>
            DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::DeterministicSchedsParetoExplorer(Environment const& env, preprocessing::SparseMultiObjectivePreprocessorResult<SparseModelType>& preprocessorResult) : model(preprocessorResult.preprocessedModel), objectives(preprocessorResult.objectives) {
                originalModelInitialState = *preprocessorResult.originalModel.getInitialStates().begin();
                lpChecker = std::make_shared<DeterministicSchedsLpChecker<SparseModelType, GeometryValueType>>(env, *model, objectives);
            }

            template <class SparseModelType, typename GeometryValueType>
            std::unique_ptr<CheckResult> DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::check(Environment const& env) {
                clean();
                initializeFacets(env);
                while (!unprocessedFacets.empty()) {
                    Facet f = std::move(unprocessedFacets.front());
                    unprocessedFacets.pop();
                    processFacet(env, f);
                }
                
                std::vector<std::vector<ModelValueType>>paretoPoints;
                paretoPoints.reserve(pointset.size());
                for (auto const& p : pointset) {
                    paretoPoints.push_back(storm::utility::vector::convertNumericVector<ModelValueType>(transformObjectiveValuesToOriginal(objectives, p.second.get())));
                }
                return std::make_unique<storm::modelchecker::ExplicitParetoCurveCheckResult<ModelValueType>>(originalModelInitialState, std::move(paretoPoints), nullptr, nullptr);
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::clean() {
                pointset = Pointset();
                unprocessedFacets = std::queue<Facet>();
                overApproximation = storm::storage::geometry::Polytope<GeometryValueType>::createUniversalPolytope();
                unachievableAreas.clear();
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::addHalfspaceToOverApproximation(Environment const& env, std::vector<GeometryValueType> const& normalVector, GeometryValueType const& offset) {
                if (env.modelchecker().multi().isPrintResultsSet()) {
                    std::cout << "## unachievable halfspace: [";
                    bool first = true;
                    for (auto const& xi : normalVector) {
                        if (first) {
                            first = false;
                        } else {
                            std::cout << ",";
                        }
                        std::cout << storm::utility::convertNumber<double>(xi);
                    }
                    std::cout << "];[" << storm::utility::convertNumber<double>(offset) << "]" << std::endl;
                }
                storm::storage::geometry::Halfspace<GeometryValueType> overApproxHalfspace(normalVector, offset);
                overApproximation = overApproximation->intersection(overApproxHalfspace);
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::addUnachievableArea(Environment const& env, Polytope const& area) {
                if (env.modelchecker().multi().isPrintResultsSet()) {
                    std::vector<std::vector<GeometryValueType>> vertices;
                    if (objectives.size() == 2) {
                        vertices = area->getVerticesInClockwiseOrder();
                    } else {
                        vertices = area->getVertices();
                    }
                    std::cout << "## unachievable polytope: ";
                    bool firstVertex = true;
                    for (auto const& v : vertices) {
                        if (firstVertex) {
                            firstVertex = false;
                        } else {
                            std::cout << ";";
                        }
                        std::cout << "[";
                        bool firstEntry = true;
                        for (auto const& vi : v) {
                            if (firstEntry) {
                                firstEntry = false;
                            } else {
                                std::cout << ",";
                            }
                            std::cout << storm::utility::convertNumber<double>(vi);
                        }
                        std::cout << "]";
                    }
                    std::cout << std::endl;
                }
                unachievableAreas.push_back(area);
            }
                
            template <class SparseModelType, typename GeometryValueType>
            typename DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::Polytope DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::negateMinObjectives(Polytope const& polytope) const {
                std::vector<GeometryValueType> zeroRow(objectives.size(), storm::utility::zero<GeometryValueType>());
                std::vector<std::vector<GeometryValueType>> transformationMatrix(objectives.size(), zeroRow);
                for (uint64_t objIndex = 0; objIndex < objectives.size(); ++objIndex) {
                    if (storm::solver::minimize(objectives[objIndex].formula->getOptimalityType())) {
                        transformationMatrix[objIndex][objIndex] = -storm::utility::one<GeometryValueType>();
                    } else {
                        transformationMatrix[objIndex][objIndex] = storm::utility::one<GeometryValueType>();
                    }
                }
                return polytope->affineTransformation(transformationMatrix, zeroRow);
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::negateMinObjectives(std::vector<GeometryValueType>& vector) const {
                for (uint64_t objIndex = 0; objIndex < this->objectives.size(); ++objIndex) {
                    if (storm::solver::minimize(objectives[objIndex].formula->getOptimalityType())) {
                        vector[objIndex] *= -storm::utility::one<ModelValueType>();
                    }
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::initializeFacets(Environment const& env) {
                for (uint64_t objIndex = 0; objIndex < objectives.size(); ++objIndex) {
                    std::vector<GeometryValueType> weightVector(objectives.size(), storm::utility::zero<ModelValueType>());
                    weightVector[objIndex] = storm::utility::one<GeometryValueType>();
                    negateMinObjectives(weightVector);
                    lpChecker->setCurrentWeightVector(weightVector);
                    auto point = lpChecker->check(env, negateMinObjectives(this->overApproximation));
                    STORM_LOG_THROW(point.is_initialized(), storm::exceptions::UnexpectedException, "Unable to find a point in the current overapproximation.");
                    negateMinObjectives(weightVector);
                    negateMinObjectives(point.get());
                    Point p(point.get());
                    p.setOnFacet();
                    // Adapt the overapproximation
                    GeometryValueType offset = storm::utility::vector::dotProduct(weightVector, p.get());
                    addHalfspaceToOverApproximation(env, weightVector, offset);
                    pointset.addPoint(env, std::move(p));
                }
                
                auto initialHalfspaces = pointset.downwardClosure()->getHalfspaces();
                for (auto& h : initialHalfspaces) {
                    Facet f(std::move(h));
                    for (auto const& p : pointset) {
                        if (f.getHalfspace().isPointOnBoundary(p.second.get())) {
                            f.addPoint(p.first, p.second);
                        }
                    }
                    STORM_LOG_ASSERT(std::count(f.getHalfspace().normalVector().begin(), f.getHalfspace().normalVector().end(), storm::utility::zero<GeometryValueType>()) + f.getNumberOfPoints() == objectives.size(), "Unexpected number of points on facet.");
                    
                    unprocessedFacets.push(std::move(f));
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            std::vector<GeometryValueType>  DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::getReferenceCoordinates() const {
                std::vector<GeometryValueType> result;
                for (auto const& obj : objectives) {
                    // TODO: use objectiveHelper here.
                    ModelValueType value = storm::solver::minimize(obj.formula->getOptimalityType()) ? -obj.upperResultBound.get() : obj.lowerResultBound.get();
                    result.push_back(storm::utility::convertNumber<GeometryValueType>(value));
                }
                return result;
            }

            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::processFacet(Environment const& env, Facet& f) {
                std::vector<GeometryValueType> weightVector = f.getHalfspace().normalVector();
                negateMinObjectives(weightVector);
                lpChecker->setCurrentWeightVector(weightVector);
                
                if (optimizeAndSplitFacet(env,f)) {
                    return;
                }
                
                GeometryValueType eps = storm::utility::convertNumber<GeometryValueType>(env.modelchecker().multi().getPrecision());
                eps += eps; // The unknown area (box) can actually have size 2*eps
                storm::storage::geometry::PolytopeTree<GeometryValueType> polytopeTree(f.getInducedSimplex(pointset, getReferenceCoordinates()));
                for (auto const& point : pointset) {
                    polytopeTree.substractDownwardClosure(point.second.get(), eps);
                    if (polytopeTree.isEmpty()) {
                        break;
                    }
                }
                if (!polytopeTree.isEmpty()) {
                    auto res = lpChecker->check(env, polytopeTree, eps);
                    for (auto const& infeasableArea : res.second) {
                        addUnachievableArea(env, negateMinObjectives(infeasableArea));
                    }
                    for (auto& achievablePoint : res.first) {
                        negateMinObjectives(achievablePoint);
                        pointset.addPoint(env, Point(std::move(achievablePoint)));
                    }
                }
            }
            
            template <class SparseModelType, typename GeometryValueType>
            bool DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::optimizeAndSplitFacet(Environment const& env, Facet& f) {
                
                // Invoke optimization and insert the explored points
                boost::optional<PointId> optPointId;
                auto currentArea = negateMinObjectives(overApproximation->intersection(f.getHalfspace().invert()));
                auto point = lpChecker->check(env, currentArea);
                if (point.is_initialized()) {
                    negateMinObjectives(point.get());
                    Point p(point.get());
                    p.setOnFacet();
                    GeometryValueType offset = storm::utility::vector::dotProduct(f.getHalfspace().normalVector(), p.get());
                    addHalfspaceToOverApproximation(env, f.getHalfspace().normalVector(), offset);
                    optPointId = pointset.addPoint(env, std::move(p));
                } else {
                    addHalfspaceToOverApproximation(env, f.getHalfspace().normalVector(), f.getHalfspace().offset());
                }
                
                // Potentially generate new facets
                if (optPointId) {
                    auto const& optPoint = pointset.getPoint(*optPointId);
                    // TODO: this check might suffer from numerical errors. Check how much this would hurt us.
                    if (f.getHalfspace().contains(optPoint.get())) {
                        // The point is contained in the halfspace which means that no more splitting is possible.
                        return false;
                    } else {
                        // Found a new Pareto optimal point -> generate new facets
                        std::vector<std::vector<GeometryValueType>> vertices;
                        vertices.push_back(optPoint.get());
                        for (auto const& pId : f.getPoints()) {
                            vertices.push_back(pointset.getPoint(pId).get());
                        }
                        auto newHalfspaceCandidates = storm::storage::geometry::Polytope<GeometryValueType>::createSelectiveDownwardClosure(vertices, storm::utility::vector::filterZero(f.getHalfspace().normalVector()))->getHalfspaces();
                        for (auto& h : newHalfspaceCandidates) {
                            if (!storm::utility::vector::hasNegativeEntry(h.normalVector())) {
                                STORM_LOG_ASSERT(h.isPointOnBoundary(optPoint.get()), "Unexpected facet found while splitting.");
                                Facet fNew(std::move(h));
                                fNew.addPoint(optPointId.get(), optPoint);
                                auto vertexIt = vertices.begin();
                                ++vertexIt;
                                for (auto const& pId : f.getPoints()) {
                                    assert(pointset.getPoint(pId).get() == *vertexIt);
                                    if (fNew.getHalfspace().isPointOnBoundary(*vertexIt)) {
                                        fNew.addPoint(pId, pointset.getPoint(pId));
                                    }
                                    ++vertexIt;
                                }
                                assert(vertexIt == vertices.end());
                                unprocessedFacets.push(std::move(fNew));
                            }
                        }
                        return true;
                    }
                } else {
                    // If the 'optimal point' was dominated by an existing point, we can not split the facet any further.
                    return false;
                }
            }

            template <class SparseModelType, typename GeometryValueType>
            void DeterministicSchedsParetoExplorer<SparseModelType, GeometryValueType>::exportPlotOfCurrentApproximation(Environment const& env) {
                /*
                STORM_LOG_ERROR_COND(objectives.size()==2, "Exporting plot requested but this is only implemented for the two-dimensional case.");
                
                auto transformedUnderApprox = transformPolytopeToOriginalModel(underApproximation);
                auto transformedOverApprox = transformPolytopeToOriginalModel(overApproximation);
                
                // Get pareto points as well as a hyperrectangle that is used to guarantee that the resulting polytopes are bounded.
                storm::storage::geometry::Hyperrectangle<GeometryValueType> boundaries(std::vector<GeometryValueType>(objectives.size(), storm::utility::zero<GeometryValueType>()), std::vector<GeometryValueType>(objectives.size(), storm::utility::zero<GeometryValueType>()));
                std::vector<std::vector<GeometryValueType>> paretoPoints;
                paretoPoints.reserve(refinementSteps.size());
                for(auto const& step : refinementSteps) {
                    paretoPoints.push_back(transformPointToOriginalModel(step.lowerBoundPoint));
                    boundaries.enlarge(paretoPoints.back());
                }
                auto underApproxVertices = transformedUnderApprox->getVertices();
                for(auto const& v : underApproxVertices) {
                    boundaries.enlarge(v);
                }
                auto overApproxVertices = transformedOverApprox->getVertices();
                for(auto const& v : overApproxVertices) {
                    boundaries.enlarge(v);
                }
                
                //Further enlarge the boundaries a little
                storm::utility::vector::scaleVectorInPlace(boundaries.lowerBounds(), GeometryValueType(15) / GeometryValueType(10));
                storm::utility::vector::scaleVectorInPlace(boundaries.upperBounds(), GeometryValueType(15) / GeometryValueType(10));
                
                auto boundariesAsPolytope = boundaries.asPolytope();
                std::vector<std::string> columnHeaders = {"x", "y"};
                
                std::vector<std::vector<double>> pointsForPlotting;
                if (env.modelchecker().multi().getPlotPathUnderApproximation()) {
                    underApproxVertices = transformedUnderApprox->intersection(boundariesAsPolytope)->getVerticesInClockwiseOrder();
                    pointsForPlotting.reserve(underApproxVertices.size());
                    for(auto const& v : underApproxVertices) {
                        pointsForPlotting.push_back(storm::utility::vector::convertNumericVector<double>(v));
                    }
                    storm::utility::exportDataToCSVFile<double, std::string>(env.modelchecker().multi().getPlotPathUnderApproximation().get(), pointsForPlotting, columnHeaders);
                }
                
                if (env.modelchecker().multi().getPlotPathOverApproximation()) {
                    pointsForPlotting.clear();
                    overApproxVertices = transformedOverApprox->intersection(boundariesAsPolytope)->getVerticesInClockwiseOrder();
                    pointsForPlotting.reserve(overApproxVertices.size());
                    for(auto const& v : overApproxVertices) {
                        pointsForPlotting.push_back(storm::utility::vector::convertNumericVector<double>(v));
                    }
                    storm::utility::exportDataToCSVFile<double, std::string>(env.modelchecker().multi().getPlotPathOverApproximation().get(), pointsForPlotting, columnHeaders);
                }
                
                if (env.modelchecker().multi().getPlotPathParetoPoints()) {
                    pointsForPlotting.clear();
                    pointsForPlotting.reserve(paretoPoints.size());
                    for(auto const& v : paretoPoints) {
                        pointsForPlotting.push_back(storm::utility::vector::convertNumericVector<double>(v));
                    }
                    storm::utility::exportDataToCSVFile<double, std::string>(env.modelchecker().multi().getPlotPathParetoPoints().get(), pointsForPlotting, columnHeaders);
                }
            };
                 */
            }

            
            template class DeterministicSchedsParetoExplorer<storm::models::sparse::Mdp<double>, storm::RationalNumber>;
            template class DeterministicSchedsParetoExplorer<storm::models::sparse::Mdp<storm::RationalNumber>, storm::RationalNumber>;
            template class DeterministicSchedsParetoExplorer<storm::models::sparse::MarkovAutomaton<double>, storm::RationalNumber>;
            template class DeterministicSchedsParetoExplorer<storm::models::sparse::MarkovAutomaton<storm::RationalNumber>, storm::RationalNumber>;
        }
    }
}