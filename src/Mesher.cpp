#include "Mesher.h"

Mesher::Mesher()
{
    inMesh_.initialize();
    outMesh_.initialize();
    info_.initialize();
}

Mesher::~Mesher() 
{
    inMesh_.finalize();
    outMesh_.finalize();
}

void Mesher::executePreMeshingProcesses(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    info_.initialize();
    setMeshInfo(nodes, elements, param);
    removeMeshNodes(nodes, elements, param);
    generateNewNodes(nodes, elements, param);
}

void Mesher::executePostMeshingProcesses(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    selectMeshElements(nodes, elements, param);
    generateNewElements(nodes, elements, param);
    buildFreeSurface(nodes, elements, param);
    setMeshInfo(nodes, elements, param);
    //deactivateSlivers(nodes, elements, param);
}

void Mesher::setMeshInfo(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    int number_of_new_nodes = nodes.size() - info_.getNumberOfNodes();
    int number_of_new_elements = elements.size() - info_.getNumberOfElements();

    if (number_of_new_nodes > 0)
        info_.setNumberOfNewNodes(number_of_new_nodes);
    else
        info_.setNumberOfNewNodes(0);

    if (number_of_new_elements > 0)
        info_.setNumberOfNewElements(number_of_new_elements);
    else
        info_.setNumberOfNewElements(0);

    info_.setNumberOfNodes(nodes.size());
    info_.setNumberOfElements(elements.size());
}


void Mesher::removeMeshNodes(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{   
    bool any_node_removed = false;
    int inside_nodes_removed = 0;
    int area_nodes_removed = 0;
    int close_wall_nodes_removed = 0;
    int sliver_wall_nodes_removed = 0;

    unsigned int dimension = param->getDimension();
    double deltat = param->getDeltat();
    double currentTime = param->getCurrentTime();
    bool firstMesh = false;
    if (currentTime < 2.0 * deltat)
    {
        firstMesh = true;
    }
    double meanMeshLength = param->getMeshLength();

    //Looking for critical elements, both slivers near to the free surface and zero area elements close to the walls, which could lead to fluid leakage
    for (Element* const& el : elements)
    {
        const std::vector<Node*>& elNodes = el->getNodes();
        int nnodes = elNodes.size();
        
        int rigidNodes = 0;
        int freeSurfaceNodes = 0;

        for (Node* const& node : elNodes)
        {
            if (node->isConstrained() || node->isInterface())
                rigidNodes++;
            freeSurfaceNodes += node->isFreeSurface();
        }

        if (rigidNodes > 1 && rigidNodes != nnodes)
        {
            double volume = el->getBaseElement()->getJacobianIntegration();
            if (dimension == 2)
            {
                double safetyCoefficient = 0.5;

                double edgeLengths[3];
                double wallLength = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                    edgeLengths[i] = elNodes[edgeVerticesId[0]]->distanceToNode(*elNodes[edgeVerticesId[1]], dimension);
                    if ((elNodes[edgeVerticesId[0]]->isConstrained() && elNodes[edgeVerticesId[1]]->isConstrained()) ||
                        (elNodes[edgeVerticesId[0]]->isInterface() && elNodes[edgeVerticesId[1]]->isInterface()) && edgeLengths[i] > wallLength)
                        wallLength = edgeLengths[i];
                }

                for (Node* const& node : elNodes)
                {
                    if (!node->isConstrained() && !node->isInterface() && !node->isToRemove())
                    { 
                        double height = 2.0 * volume / wallLength;

                        if (node->isFreeSurface())
                        {
                            const std::vector<Node*>& neighborNodes = node->getNeighborNodes();
                            int neighborRigid = 0;
                            int neighborFreeSurface = 0;
                            for (Node* const& neighborNode : neighborNodes)
                            {
                                if (neighborNode->isConstrained() || neighborNode->isInterface())
                                    neighborRigid++;
                                else if (neighborNode->isFreeSurface())
                                    neighborFreeSurface++;
                            }
                            if ((neighborRigid + neighborFreeSurface) == neighborNodes.size())
                                safetyCoefficient = 0.25;
                        }

                        if (height < (0.5 * safetyCoefficient * wallLength))
                        {
                            node->setToRemove(true);
                            any_node_removed = true;
                            close_wall_nodes_removed++;
                            break;
                        }
                    }
                }

            }
            else if (dimension == 3)
            {
                double safetyCoefficient = 0.6;
                double criticalVolume = 0.1 * param->getModelVolume() / (double)elements.size();

                if (volume < criticalVolume)
                {
                    for (Node* const& node : elNodes)
                    {
                        if (!node->isConstrained() && !node->isInterface() && !node->isToRemove())
                        {
                            node->setToRemove(true);
                            any_node_removed = true;
                            area_nodes_removed++;
                        }
                    }
                }

                else if (rigidNodes == 3)
                {
                    int rigidNodesId[3];
                    int notRigidNodeId;
                    int aux = -1;
                    for (int i = 0; i < nnodes; i++)
                    {
                        if (elNodes[i]->isConstrained() || elNodes[i]->isInterface())
                            rigidNodesId[++aux]  = i;
                        else
                            notRigidNodeId = i;
                    }

                    if (!elNodes[notRigidNodeId]->isToRemove())
                    {
                        double a1; //slope x for the plane composed by rigid nodes only
                        double b1; //slope y for the plane composed by rigid nodes only
                        double c1; //slope z for the plane composed by rigid nodes only
                        a1 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) - 
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue());
                        b1 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) - 
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue());
                        c1 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) -
                             (elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue());
                        double a2; //slope x for the plane composed by rigid nodes 0, 1 and the non rigid node
                        double b2; //slope y for the plane composed by rigid nodes 0, 1 and the non rigid node
                        double c2; //slope z for the plane composed by rigid nodes 0, 1 and the non rigid node
                        a2 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue());
                        b2 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue());
                        c2 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue());
                        double a3; //slope x for the plane composed by rigid nodes 1, 2 and the non rigid node
                        double b3; //slope y for the plane composed by rigid nodes 1, 2 and the non rigid node
                        double c3; //slope z for the plane composed by rigid nodes 1, 2 and the non rigid node
                        a3 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue());
                        b3 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue());
                        c3 = (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[rigidNodesId[1]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue());
                        double a4; //slope x for the plane composed by rigid nodes 0, 2 and the non rigid node
                        double b4; //slope y for the plane composed by rigid nodes 0, 2 and the non rigid node
                        double c4; //slope z for the plane composed by rigid nodes 0, 2 and the non rigid node
                        a4 = (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                             (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue());
                        b4 = (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(2)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(2)->getCurrentValue()) *
                             (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue());
                        c4 = (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(1)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue()) -
                             (elNodes[notRigidNodeId]->getDegreeOfFreedom(0)->getCurrentValue()  - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(0)->getCurrentValue()) *
                             (elNodes[rigidNodesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[rigidNodesId[2]]->getDegreeOfFreedom(1)->getCurrentValue());

                        //angle between the rigid plane and the other plans. If the angle is small, the particle can pass through the wall
                        double factor = sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2));
                        double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (factor * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                        double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (factor * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
                        double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (factor * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));

                        if (fabs(cosAngle12) > 0.995 || fabs(cosAngle13) > 0.995 || fabs(cosAngle14) > 0.995)
                        {
                            elNodes[notRigidNodeId]->setToRemove(true);
                            any_node_removed = true;
                            sliver_wall_nodes_removed++;
                        }
                    }
                }

                double edgeLengths[6];
                double wallLength = 0.0;
                for (int i = 0; i < 6; i++)
                {
                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                    edgeLengths[i] = elNodes[edgeVerticesId[0]]->distanceToNode(*elNodes[edgeVerticesId[1]], dimension);
                    if (((elNodes[edgeVerticesId[0]]->isConstrained() && elNodes[edgeVerticesId[1]]->isConstrained()) ||
                        (elNodes[edgeVerticesId[0]]->isInterface() && elNodes[edgeVerticesId[1]]->isInterface())) && wallLength == 0.0)
                        wallLength = edgeLengths[i];
                }

                for (Node* const& node : elNodes)
                {
                    if (!node->isConstrained() && !node->isInterface() && !node->isIsolated() && !node->isToRemove() && node->isFreeSurface())
                    {
                        const std::vector<Node*>& neighbor_nodes = node->getNeighborNodes();
                        int neighborRigid = 0;
                        int neighborFreeSurface = 0;
                        for (Node* const& neighbor_node : neighbor_nodes)
                        {
                            if (neighbor_node->isConstrained() || neighbor_node->isInterface())
                                neighborRigid++;
                            else if (neighbor_node->isFreeSurface())
                                neighborFreeSurface++;
                        }
                        if ((neighborRigid + neighborFreeSurface) == neighbor_nodes.size() && neighborRigid > 0)
                            safetyCoefficient = 0.25;
                    }
                }

                for (int i = 0; i < 6; i++)
                {
                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                    bool firstNodeRigid = false, secondNodeRigid = false;
                    if (elNodes[edgeVerticesId[0]]->isConstrained() || elNodes[edgeVerticesId[0]]->isInterface())
                        firstNodeRigid = true;
                    if (elNodes[edgeVerticesId[1]]->isConstrained() || elNodes[edgeVerticesId[1]]->isInterface())
                        secondNodeRigid = true;
                    if (((firstNodeRigid && !secondNodeRigid) || (!firstNodeRigid && secondNodeRigid)) &&
                        !elNodes[edgeVerticesId[0]]->isToRemove() && !elNodes[edgeVerticesId[1]]->isToRemove() &&
                        edgeLengths[i] < wallLength * safetyCoefficient)
                    {
                        if (!firstNodeRigid && !elNodes[edgeVerticesId[0]]->isToRemove())
                        {
                            elNodes[edgeVerticesId[0]]->setToRemove(true);
                            any_node_removed = true;
                            close_wall_nodes_removed++;
                        }
                        else if (!secondNodeRigid && !elNodes[edgeVerticesId[1]]->isToRemove())
                        {
                            elNodes[edgeVerticesId[1]]->setToRemove(true);
                            any_node_removed = true;
                            close_wall_nodes_removed++;
                        }
                    }
                }
            }
        }
    }

    int threshold_distance_boundary = 0.6 * meanMeshLength;

    //Check if there are nodes too close to each other
    for (Node* const& node : nodes)
    {
        if (!node->isNewEntity() && !node->isIsolated() && !node->isToRemove())
        {
            const unsigned int ndofs = node->getNumberOfDegreesOfFreedom();
            double radius = 0.6 * meanMeshLength;
            unsigned int neighborRemovedNodes = 0;
            unsigned int freeSurfaceNeighborNodes = 0;
            const std::vector<Node*>& neighborNodes = node->getNeighborNodes();
            auto it_begin = neighborNodes.begin() + 1; //the first neighbor node is always the node itself
            auto it_end = neighborNodes.end();

            if (node->isFreeSurface()) //if a node belongs to the free surface, it must be more difficult to remove it. Otherwise, a lot of volume is lost
            {
                radius *= 0.75;
                unsigned int neighborRigids = 0;
                
                for (auto& nnode = it_begin; nnode != it_end; nnode++)
                {
                    if ((*nnode)->isConstrained() || (*nnode)->isInterface())
                        neighborRigids++;
                    if ((*nnode)->isToRemove())
                        neighborRemovedNodes++;
                }
                if (neighborRigids == (neighborNodes.size() - 1))
                    radius *= 0.3;

            }
            else
            {
                for (auto& nnode = it_begin; nnode != it_end; nnode++)
                {
                    if ((*nnode)->isFreeSurface())
                        freeSurfaceNeighborNodes++;
                    if ((*nnode)->isToRemove())
                        neighborRemovedNodes++;
                }
            }

            if (freeSurfaceNeighborNodes > 1)
                radius = 0.5 * meanMeshLength;
            
            int nFoundNodes; std::vector<Node*> foundNodes; std::vector<double> distances;
            node->searchNodesInRadius(radius, dimension, nFoundNodes, foundNodes, distances);

            if (nFoundNodes > 1 && neighborRemovedNodes == 0)
            {
                if (!node->isConstrained() && !node->isInterface())
                {
                    if (!node->isFreeSurface() && freeSurfaceNeighborNodes == dimension) // if a node is close to the free surface, we move it back rather than erasing it
                    {
                        for (int dof = 0; dof < ndofs; dof++)
                        {
                            double val = 0.0;
                            double firstDerivative = 0.0;
                            double secondDerivative = 0.0;
                            for (Node* const& nnode : neighborNodes)
                            {
                                val += nnode->getDegreeOfFreedom(dof)->getCurrentValue();
                                firstDerivative += nnode->getDegreeOfFreedom(dof)->getCurrentFirstTimeDerivative();
                                secondDerivative += nnode->getDegreeOfFreedom(dof)->getCurrentSecondTimeDerivative();
                            }
                            node->getDegreeOfFreedom(dof)->setCurrentValue(val/neighborNodes.size());
                            node->getDegreeOfFreedom(dof)->setCurrentFirstTimeDerivative(firstDerivative/neighborNodes.size());
                            node->getDegreeOfFreedom(dof)->setCurrentSecondTimeDerivative(secondDerivative/neighborNodes.size());
                        }
                    }
                    else
                    {
                        node->setToRemove(true);
                        any_node_removed = true;
                        inside_nodes_removed++;
                    }
                }
            }
        }
    }

    info_.removedNodes_ = inside_nodes_removed + area_nodes_removed + close_wall_nodes_removed + sliver_wall_nodes_removed;

    // We cant delete the nodes marked to be erased here,
    // otherwise the pointers inside the elements would be invalidated when accessing it in generateNewNodes method.
}

void Mesher::generateNewNodes(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    int dimension = param->getDimension();
    double deltat = param->getDeltat();
    double currentTime = param->getCurrentTime();

    if (currentTime < 2.0 * deltat)
    {
        //info_.removedNodes_ = 0;
        info_.numberOfInitialNodes_ = info_.numberOfNodes_;
    }

    int elementsToRefine = info_.removedNodes_; //We try to refine the same number of elements as the nodes removed
    int extraNodes = info_.numberOfNodes_ - info_.numberOfInitialNodes_;
    int toleratedExtraNodes = int(0.05 * info_.numberOfInitialNodes_);

    if (elementsToRefine > 0)
    {
        std::vector<double> biggestVolumes(elementsToRefine, -1.0);
        std::vector<std::vector<Node*>> nodesToInterpolate(elementsToRefine, std::vector<Node*>(2, nullptr));
        std::vector<std::vector<double>> newPositions(elementsToRefine, std::vector<double>(dimension));
        int countNodes = 0; //the number of new Nodes that will in fact be created

        // searching for elementsToRefine large edges to be refined
        for (Element* const& el : elements)
        {
            const std::vector<Node*>& elNodes = el->getNodes();

            int numRigid = 0;
            int numBoundary = 0;
            int numFreeSurface = 0;
            bool hasRemovedNode = false;

            for (Node* const& node : elNodes)
            {
                numRigid += node->isConstrained();
                numBoundary += node->isBoundary();
                numFreeSurface += node->isFreeSurface();
                if (node->isToRemove()) hasRemovedNode = true;
            }
            
            if (dimension == 2)
            {
                double edgeLengthTreshold = 1.4 * param->getMeshLength();
                double safetyCoefficient = 1.5;
                double volume = el->getBaseElement()->getJacobianIntegration();

                double edgeLengths[3];
                double wallLength = 0.0;
                for (int i = 0; i < 3; i++)
                {
                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                    edgeLengths[i] = elNodes[edgeVerticesId[0]]->distanceToNode(*elNodes[edgeVerticesId[1]], dimension);
                    if (elNodes[edgeVerticesId[0]]->isConstrained() && elNodes[edgeVerticesId[1]]->isConstrained() && edgeLengths[i] > wallLength)
                        wallLength = edgeLengths[i];
                }

                bool dangerousElementToRefine = false;
                if (numRigid > 1)
                {
                    for (unsigned int i = 0; i < 3; i++)
                    {
                        const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                        if (edgeLengths[i] < wallLength * safetyCoefficient & (elNodes[edgeVerticesId[0]]->isConstrained() || elNodes[edgeVerticesId[1]]->isConstrained()))
                            edgeLengths[i] = 0.0;
                        else if ((elNodes[edgeVerticesId[0]]->isFreeSurface() || elNodes[edgeVerticesId[0]]->isConstrained()) &&
                                 (elNodes[edgeVerticesId[1]]->isFreeSurface() || elNodes[edgeVerticesId[1]]->isConstrained()))
                            edgeLengths[i] = 0.0;
                    }
                }
                if (edgeLengths[0] == 0.0 && edgeLengths[1] == 0.0 && edgeLengths[2] == 0.0)
                    dangerousElementToRefine = true;
                
                if (!dangerousElementToRefine && !hasRemovedNode)
                {
                    auto it = std::max_element(edgeLengths, edgeLengths+3);
                    int largestEdge = std::distance(edgeLengths, it);
                    double largestEdgeLength = *it;

                    if (countNodes < elementsToRefine && largestEdgeLength > edgeLengthTreshold)
                    {
                        const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(largestEdge);
                        double newPosition[2];
                        newPosition[0] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() +
                                          elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(0)->getCurrentValue()) * 0.5;
                        newPosition[1] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() +
                                          elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(1)->getCurrentValue()) * 0.5;
                        
                        bool suitableElement = true;
                        for (int i = 0; i < countNodes; i++)
                        {
                            double diff_x = fabs(newPositions[i][0] - newPosition[0]) - param->getMeshLength() * 0.5;
                            double diff_y = fabs(newPositions[i][1] - newPosition[1]) - param->getMeshLength() * 0.5;
                            if (diff_x < 0.0 && diff_y < 0.0) // the considered new node is in the same region of a previously generated node
                                suitableElement = false;
                        }

                        if (suitableElement)
                        {
                            nodesToInterpolate[countNodes][0] = elNodes[edgeVerticesId[0]];
                            nodesToInterpolate[countNodes][1] = elNodes[edgeVerticesId[1]];
                            biggestVolumes[countNodes] = volume;
                            newPositions[countNodes][0] = newPosition[0];
                            newPositions[countNodes][1] = newPosition[1];
                            countNodes++;
                        }
                    }
                    else if (numFreeSurface < 3)
                    {
                        double penalization = 1.0;
                        if (numRigid > 1)
                            penalization = 0.8;
                        else if (numRigid > 0 && numFreeSurface > 0)
                            penalization = 0.0;
                        else if (numFreeSurface > 0)
                            penalization = 0.875;
                        volume *= penalization;

                        for (int i = 0; i < elementsToRefine; i++)
                        {
                            if (volume > biggestVolumes[i])
                            {
                                bool suitableElement = true;
                                if (largestEdgeLength > edgeLengthTreshold)
                                {
                                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(largestEdge);
                                    double newPosition[2];
                                    newPosition[0] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() +
                                                    elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(0)->getCurrentValue()) * 0.5;
                                    newPosition[1] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() +
                                                    elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(1)->getCurrentValue()) * 0.5;
                                    
                                    bool suitableElement = true;
                                    for (int j = 0; j < countNodes; j++)
                                    {
                                        double diff_x = fabs(newPositions[j][0] - newPosition[0]) - param->getMeshLength() * 0.5;
                                        double diff_y = fabs(newPositions[j][1] - newPosition[1]) - param->getMeshLength() * 0.5;
                                        if (diff_x < 0.0 && diff_y < 0.0) // the considered new node is in the same region of a previously generated node
                                            suitableElement = false;
                                    }

                                    if (suitableElement)
                                    {
                                        nodesToInterpolate[i][0] = elNodes[edgeVerticesId[0]];
                                        nodesToInterpolate[i][1] = elNodes[edgeVerticesId[1]];
                                        biggestVolumes[i] = volume;
                                        newPositions[i][0] = newPosition[0];
                                        newPositions[i][1] = newPosition[1];
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
                
            }
            else if (dimension == 3)
            {
                double edgeLengthTreshold = 1.25 * param->getMeshLength();
                double safetyCoefficient = 1.6;
                double volume = el->getBaseElement()->getJacobianIntegration();

                double edgeLengths[6];
                double wallLength = 0.0;
                for (int i = 0; i < 6; i++)
                {
                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                    edgeLengths[i] = elNodes[edgeVerticesId[0]]->distanceToNode(*elNodes[edgeVerticesId[1]], dimension);
                    if (elNodes[edgeVerticesId[0]]->isConstrained() && elNodes[edgeVerticesId[1]]->isConstrained() && edgeLengths[i] > wallLength)
                        wallLength = edgeLengths[i];
                }

                bool dangerousElementToRefine = false;
                if (numRigid == 1)
                {
                    if (elNodes[0]->isConstrained())
                    {
                        edgeLengths[0] = 0.0;
                        edgeLengths[3] = 0.0;
                        edgeLengths[5] = 0.0;
                    }
                    else if (elNodes[1]->isConstrained())
                    {
                        edgeLengths[0] = 0.0;
                        edgeLengths[1] = 0.0;
                        edgeLengths[4] = 0.0;
                    }
                    else if (elNodes[2]->isConstrained())
                    {
                        edgeLengths[1] = 0.0;
                        edgeLengths[2] = 0.0;
                        edgeLengths[5] = 0.0;
                    }
                    else if (elNodes[3]->isConstrained())
                    {
                        edgeLengths[2] = 0.0;
                        edgeLengths[3] = 0.0;
                        edgeLengths[4] = 0.0;
                    }
                }
                else if (numRigid == 2)
                {
                    for (unsigned int i = 0; i < 6; i++)
                    {
                        const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(i);
                        if (edgeLengths[i] < wallLength * safetyCoefficient && (elNodes[edgeVerticesId[0]]->isConstrained() || 
                            elNodes[edgeVerticesId[1]]->isConstrained()))
                            edgeLengths[i] = 0.0;
                        else if ((elNodes[edgeVerticesId[0]]->isFreeSurface() || elNodes[edgeVerticesId[0]]->isConstrained()) &&
                                 (elNodes[edgeVerticesId[1]]->isFreeSurface() || elNodes[edgeVerticesId[1]]->isConstrained()))
                            edgeLengths[i] = 0.0;
                    }
                }

                if ((edgeLengths[0] == 0.0 && edgeLengths[1] == 0.0 && edgeLengths[2] == 0.0 &&
                    edgeLengths[3] == 0.0 && edgeLengths[4] == 0.0 && edgeLengths[5] == 0.0) || numRigid > 2)
                    dangerousElementToRefine = true;

                if (!dangerousElementToRefine && !hasRemovedNode)
                {
                    auto it = std::max_element(edgeLengths, edgeLengths+6);
                    int largestEdge = std::distance(edgeLengths, it);
                    double largestEdgeLength = *it;

                    if (countNodes < elementsToRefine && largestEdgeLength > edgeLengthTreshold)
                    {
                        const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(largestEdge);
                        double newPosition[3];
                        newPosition[0] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() +
                                          elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(0)->getCurrentValue()) * 0.5;
                        newPosition[1] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() +
                                          elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(1)->getCurrentValue()) * 0.5;
                        newPosition[2] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(2)->getCurrentValue() +
                                          elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(2)->getCurrentValue()) * 0.5;
                        
                        bool suitableElement = true;
                        for (int i = 0; i < countNodes; i++)
                        {
                            double diff_x = fabs(newPositions[i][0] - newPosition[0]) - param->getMeshLength() * 0.5;
                            double diff_y = fabs(newPositions[i][1] - newPosition[1]) - param->getMeshLength() * 0.5;
                            double diff_z = fabs(newPositions[i][2] - newPosition[2]) - param->getMeshLength() * 0.5;
                            if (diff_x < 0.0 && diff_y < 0.0 && diff_z < 0.0) // the considered new node is in the same region of a previously generated node
                                suitableElement = false;
                        }

                        if (suitableElement)
                        {
                            nodesToInterpolate[countNodes][0] = elNodes[edgeVerticesId[0]];
                            nodesToInterpolate[countNodes][1] = elNodes[edgeVerticesId[1]];
                            biggestVolumes[countNodes] = volume;
                            newPositions[countNodes][0] = newPosition[0];
                            newPositions[countNodes][1] = newPosition[1];
                            newPositions[countNodes][2] = newPosition[2];
                            countNodes++;
                        }
                    }
                    else if (numFreeSurface < 4)
                    {
                        double penalization = 1.0;
                        if (numRigid > 2)
                            penalization = 0.7;
                        else if (numRigid > 0 && numFreeSurface > 0)
                            penalization = 0.0;
                        else if (numFreeSurface > 0)
                            penalization = 0.95;
                        volume *= penalization;

                        for (int i = 0; i < elementsToRefine; i++)
                        {
                            if (volume > biggestVolumes[i])
                            {
                                bool suitableElement = true;
                                if (largestEdgeLength > edgeLengthTreshold)
                                {
                                    const std::vector<int>& edgeVerticesId = el->getParametricElement()->getEdgeVertices(largestEdge);
                                    double newPosition[3];
                                    newPosition[0] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(0)->getCurrentValue() +
                                                    elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(0)->getCurrentValue()) * 0.5;
                                    newPosition[1] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(1)->getCurrentValue() +
                                                    elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(1)->getCurrentValue()) * 0.5;
                                    newPosition[2] = (elNodes[edgeVerticesId[0]]->getDegreeOfFreedom(2)->getCurrentValue() +
                                                    elNodes[edgeVerticesId[1]]->getDegreeOfFreedom(2)->getCurrentValue()) * 0.5;
                                    
                                    bool suitableElement = true;
                                    for (int j = 0; j < countNodes; j++)
                                    {
                                        double diff_x = fabs(newPositions[j][0] - newPosition[0]) - param->getMeshLength() * 0.5;
                                        double diff_y = fabs(newPositions[j][1] - newPosition[1]) - param->getMeshLength() * 0.5;
                                        double diff_z = fabs(newPositions[j][2] - newPosition[2]) - param->getMeshLength() * 0.5;
                                        if (diff_x < 0.0 && diff_y < 0.0 && diff_z < 0.0) // the considered new node is in the same region of a previously generated node
                                            suitableElement = false;
                                    }

                                    if (suitableElement)
                                    {
                                        nodesToInterpolate[i][0] = elNodes[edgeVerticesId[0]];
                                        nodesToInterpolate[i][1] = elNodes[edgeVerticesId[1]];
                                        biggestVolumes[i] = volume;
                                        newPositions[i][0] = newPosition[0];
                                        newPositions[i][1] = newPosition[1];
                                        newPositions[i][2] = newPosition[2];
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }
        }

        //Remove and delete the nodes marked with toRemove tag (assigned in removeMeshNodes)
        std::vector<Node*> temp_nodes; temp_nodes.reserve(nodes.size());
        temp_nodes.swap(nodes);

        int node_id = -1;
        for (Node*& node : temp_nodes)
        {
            if (!node->isToRemove())
            {
                node->setIndex(++node_id);
                nodes.push_back(node);
            }
            else
            {
                delete node; //The pointers pointing to node will be invalidated inside other classes like Element
            }
        }
        temp_nodes.clear();

        info_.removedNodes_ -= countNodes;

        node_id = nodes.size() - 1;
        //Creating new countNodes nodes 
        for (int i = 0; i < countNodes; i++)
        {
           int ndofs = nodesToInterpolate[i][0]->getNumberOfDegreesOfFreedom();
           Material* mat = (nodesToInterpolate[i][0]->getMaterial()) ?
                           nodesToInterpolate[i][0]->getMaterial() : nodesToInterpolate[i][1]->getMaterial();
           std::vector<DegreeOfFreedom*> dofs; dofs.reserve(ndofs);
           for (int j = 0; j < dimension; j++)
           {
                double value = (nodesToInterpolate[i][0]->getDegreeOfFreedom(j)->getCurrentValue() +
                               nodesToInterpolate[i][1]->getDegreeOfFreedom(j)->getCurrentValue()) * 0.5;
                DegreeOfFreedom* dof = new DegreeOfFreedom(DOFType::POSITION, value);
                value = (nodesToInterpolate[i][0]->getDegreeOfFreedom(j)->getCurrentFirstTimeDerivative() +
                        nodesToInterpolate[i][1]->getDegreeOfFreedom(j)->getCurrentFirstTimeDerivative()) * 0.5;
                dof->setCurrentFirstTimeDerivative(value);
                value = (nodesToInterpolate[i][0]->getDegreeOfFreedom(j)->getCurrentSecondTimeDerivative() +
                        nodesToInterpolate[i][1]->getDegreeOfFreedom(j)->getCurrentSecondTimeDerivative()) * 0.5;
                dof->setCurrentSecondTimeDerivative(value);
                dofs.push_back(dof);
           }
           for (int j = dimension; j < ndofs; j++)
           {
                double value = (nodesToInterpolate[i][0]->getDegreeOfFreedom(j)->getCurrentValue() +
                               nodesToInterpolate[i][1]->getDegreeOfFreedom(j)->getCurrentValue()) * 0.5;
                DegreeOfFreedom* dof = new DegreeOfFreedom(DOFType::PRESSURE, value);
                dofs.push_back(dof);
           }
           Node* n = new Node(++node_id, dofs);
           n->setMaterial(mat);
           nodes.push_back(n);
        }
    }
}

void Mesher::selectMeshElements(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    // IMPORTANT!! When this function is called, neither the free surface and the isolated particles were identified.
    //             Thus, the tags isFreeSurface, isIsolated refer to the old mesh.

    int dimension = param->getDimension();

    int* outElementList = outMesh_.getElementList();
    int outNumberOfElements = outMesh_.getNumberOfElements();

    preservedElements_.clear();
    preservedElements_.resize(outNumberOfElements, -1);
    preservedElements_.shrink_to_fit();

    double deltat = param->getDeltat();
    double currentTime = param->getCurrentTime();
    bool firstMesh = false;
    if (currentTime < 2.0 * deltat)
    {
        firstMesh = true;
    }

    double modelVolume = param->getModelVolume();
    double criticalVolume = 0.05 * modelVolume / (double)elements.size();

    double meanMeshLength = param->getMeshLength();
    unsigned int numberOfElementNodes = elements[0]->getNodes().size(); //assuming that the element type does not change with remesh
    int numberOfSlivers = 0;
    int numberOfAcceptedElements = -1;

    for (int el = 0; el < outNumberOfElements; el++)
    {
        double alpha = param->getAlpha();
        
        unsigned int numFreeSurface = 0;
        unsigned int numRigid = 0;
        unsigned int numIsolated = 0;
        unsigned int numPrevFreeSurface = 0;
        unsigned int numPrevIsolated = 0;
        unsigned int sumIsolatedFreeSurface = 0;
        unsigned int sumPrevIsolatedFreeSurface = 0;
        unsigned int numIsolatedInElement = 0;
        unsigned int checkedNodes = 0;
        std::vector<double> normVelocity(numberOfElementNodes);
        std::vector<std::vector<double>> nodesVelocities(numberOfElementNodes, std::vector<double>(dimension));
        
        std::vector<Node*> elNodes(numberOfElementNodes);
        for (unsigned int i = 0; i < numberOfElementNodes; i++)
        {
            elNodes[i] = nodes[outElementList[el * numberOfElementNodes + i]];

            if (!elNodes[i]->isConstrained() && !elNodes[i]->isInterface())
            {
                const std::vector<Node*>& neighbors = elNodes[i]->getNeighborNodes();
                unsigned int nneighbors = neighbors.size() - 1;
                if (nneighbors < dimension + 1)
                    numIsolatedInElement++;
            }

            if (elNodes[i]->isConstrained() || elNodes[i]->isInterface())
            {
                numRigid++;
            }
            if (!elNodes[i]->isConstrained() && elNodes[i]->isFreeSurface())
            {
                numFreeSurface++;
                double velNorm = 0.0;
                for (int j = 0; j < dimension; j++)
                {
                    nodesVelocities[i][j] = elNodes[i]->getDegreeOfFreedom(j)->getCurrentFirstTimeDerivative();
                    velNorm += nodesVelocities[i][j] * nodesVelocities[i][j];
                }
                normVelocity[i] = sqrt(velNorm);
                checkedNodes++;
            }
            else if (elNodes[i]->isIsolated())
            {
                numIsolated++;
                double velNorm = 0.0;
                for (int j = 0; j < dimension; j++)
                {
                    nodesVelocities[i][j] = elNodes[i]->getDegreeOfFreedom(j)->getCurrentFirstTimeDerivative();
                    velNorm += nodesVelocities[i][j] * nodesVelocities[i][j];
                }
                normVelocity[i] = sqrt(velNorm);
                checkedNodes++;
            }
            if (elNodes[i]->wasFreeSurface())
                numPrevFreeSurface++;
            if (elNodes[i]->wasIsolated())
                numPrevIsolated++;
        }

        sumIsolatedFreeSurface = numIsolated + numFreeSurface;
        sumPrevIsolatedFreeSurface = numPrevFreeSurface + numPrevIsolated;

        if (dimension == 2)
        {
            if (numRigid == 0 && sumIsolatedFreeSurface == 0 && sumPrevIsolatedFreeSurface == 0)
            {
                alpha *= 1.5;
            }
            else if (sumIsolatedFreeSurface == 0 && sumPrevIsolatedFreeSurface == 0)
            {
                alpha *= 1.25;
            }
            else if (numIsolated == 0 && numPrevIsolated == 0 && numFreeSurface < numberOfElementNodes && numPrevFreeSurface < numberOfElementNodes)
            {
                alpha *= 1.125;
            }
            else
            {
                alpha *= 0.975;
            }     
        }
        else if (dimension == 3)
        {
            if (numRigid == 0 && sumIsolatedFreeSurface == 0 && sumPrevIsolatedFreeSurface == 0)
            {
                alpha *= 1.5;
            }
            else if (sumIsolatedFreeSurface == 0 && sumPrevIsolatedFreeSurface == 0)
            {
                alpha *= 1.25;
            }
            else if (numIsolated == 0 && numPrevIsolated == 0 && numFreeSurface < numberOfElementNodes && numPrevFreeSurface < numberOfElementNodes)
            {
                alpha *= 1.05;
            }
            else
            {
                alpha *= 0.95;
            } 
        }
        if (firstMesh)
        {
            alpha *= 1.15;
        }

        bool accepted = alphaShape(elNodes, alpha, meanMeshLength);

        if (numRigid == numberOfElementNodes)
        {
            accepted = false;
        }

        if (accepted && !firstMesh && (numFreeSurface == numberOfElementNodes || sumIsolatedFreeSurface == numberOfElementNodes || sumPrevIsolatedFreeSurface == numberOfElementNodes))
        {
            if (dimension == 2)
            {
                if (sumIsolatedFreeSurface == numberOfElementNodes && numRigid == 0 && checkedNodes == numberOfElementNodes)
                {
                    const double maxValue = 1.5;
                    const double minValue = 1.0 / maxValue;
                    if (normVelocity[0] / normVelocity[1] > maxValue || normVelocity[0] / normVelocity[1] < minValue ||
                        normVelocity[0] / normVelocity[2] > maxValue || normVelocity[0] / normVelocity[2] < minValue ||
                        normVelocity[1] / normVelocity[2] > maxValue || normVelocity[1] / normVelocity[2] < minValue)
                    {
                        accepted = false;
                    }
                    else
                    {
                        double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1]) /
                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2)) *
                                             sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2)));
                        double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1]) /
                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2)) *
                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2)));
                        double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1]) /
                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2)) *
                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2)));

                        if (fabs(cosAngle01) < 0.95 || fabs(cosAngle02) < 0.95 || fabs(cosAngle12) < 0.95)
                        {
                            accepted = false;
                        }
                    }
                }
                BaseSurfaceElement* triangle = new BaseSurfaceElement(0, ParametricSurfaceElement::T3, elNodes);
                double area = triangle->getJacobianIntegration();
                if (area < criticalVolume)
                {
                    accepted = false;
                }
                delete triangle;
            }
            else if (dimension == 3)
            {
                if ((sumIsolatedFreeSurface == numberOfElementNodes || numPrevIsolated == numberOfElementNodes || numPrevFreeSurface == numberOfElementNodes) 
                     && numRigid == 0 && numIsolatedInElement > 0 && checkedNodes == numberOfElementNodes)
                {
                    const double maxValue = 2.5;
                    const double minValue = 1.0 / maxValue;
                    if (normVelocity[0] / normVelocity[1] < minValue || normVelocity[0] / normVelocity[2] < minValue || normVelocity[0] / normVelocity[3] < minValue ||
                        normVelocity[0] / normVelocity[1] > maxValue || normVelocity[0] / normVelocity[2] > maxValue || normVelocity[0] / normVelocity[3] > maxValue ||
                        normVelocity[1] / normVelocity[2] < minValue || normVelocity[1] / normVelocity[3] < minValue ||
                        normVelocity[1] / normVelocity[2] > maxValue || normVelocity[1] / normVelocity[3] > maxValue ||
                        normVelocity[2] / normVelocity[3] < minValue ||
                        normVelocity[2] / normVelocity[3] > maxValue)
                    {
                        accepted = false;
                    }
                    else
                    {
                        double cosAngle01 = (nodesVelocities[0][0] * nodesVelocities[1][0] + nodesVelocities[0][1] * nodesVelocities[1][1] + nodesVelocities[0][1] * nodesVelocities[1][2]) /
                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                             sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)));
                        double cosAngle02 = (nodesVelocities[0][0] * nodesVelocities[2][0] + nodesVelocities[0][1] * nodesVelocities[2][1] + nodesVelocities[0][1] * nodesVelocities[2][2]) /
                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)));
                        double cosAngle03 = (nodesVelocities[0][0] * nodesVelocities[3][0] + nodesVelocities[0][1] * nodesVelocities[3][1] + nodesVelocities[0][1] * nodesVelocities[3][2]) /
                                            (sqrt(pow(nodesVelocities[0][0], 2) + pow(nodesVelocities[0][1], 2) + pow(nodesVelocities[0][2], 2)) *
                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));
                        double cosAngle12 = (nodesVelocities[1][0] * nodesVelocities[2][0] + nodesVelocities[1][1] * nodesVelocities[2][1] + nodesVelocities[1][1] * nodesVelocities[2][2]) /
                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)) *
                                             sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)));
                        double cosAngle13 = (nodesVelocities[1][0] * nodesVelocities[3][0] + nodesVelocities[1][1] * nodesVelocities[3][1] + nodesVelocities[1][1] * nodesVelocities[3][2]) /
                                            (sqrt(pow(nodesVelocities[1][0], 2) + pow(nodesVelocities[1][1], 2) + pow(nodesVelocities[1][2], 2)) *
                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));
                        double cosAngle23 = (nodesVelocities[2][0] * nodesVelocities[3][0] + nodesVelocities[2][1] * nodesVelocities[3][1] + nodesVelocities[2][1] * nodesVelocities[3][2]) /
                                            (sqrt(pow(nodesVelocities[2][0], 2) + pow(nodesVelocities[2][1], 2) + pow(nodesVelocities[2][2], 2)) *
                                             sqrt(pow(nodesVelocities[3][0], 2) + pow(nodesVelocities[3][1], 2) + pow(nodesVelocities[3][2], 2)));

                        if (fabs(cosAngle01) < 0.85 || fabs(cosAngle02) < 0.85 || fabs(cosAngle03) < 0.85 || fabs(cosAngle12) < 0.85 || fabs(cosAngle13) < 0.85 || fabs(cosAngle23) < 0.85)
                        {
                            accepted = false;
                        }
                    }
                }
            }
        }

        if (accepted && dimension == 3 & numRigid < 3 && (numPrevIsolated == 4 || numPrevFreeSurface == 4 || sumIsolatedFreeSurface == 4 || numFreeSurface == 4 || numIsolated == 4 || (numRigid == 2 && numIsolatedInElement > 1)))
        {
            BaseVolumeElement* tetrahedron = new BaseVolumeElement(0, ParametricVolumeElement::TET4, elNodes);
            double volume = tetrahedron->getJacobianIntegration();
            delete tetrahedron;

            if (volume < criticalVolume)
            {
                accepted = false;
                numberOfSlivers++;
            }
            else
            {
                double a1; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                double b1; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                double c1; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                a1 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) *
                     (elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) - 
                     (elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue());
                b1 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) - 
                     (elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue());
                c1 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue());
                double a2; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                double b2; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                double c2; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                a2 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) *
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue());
                b2 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue());
                c2 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue());
                double a3; //slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                double b3; //slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                double c3; //slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                a3 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue());
                b3 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue());
                c3 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue());
                double a4; //slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                double b4; //slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                double c4; //slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                a4 = (elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue());
                b4 = (elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue());
                c4 = (elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue());

                double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
                double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
                double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

                // if ((fabs(cosAngle12) > 0.99 || fabs(cosAngle13) > 0.99 || fabs(cosAngle14) > 0.99 || fabs(cosAngle23) > 0.99 ||
                //      fabs(cosAngle24) > 0.99 || fabs(cosAngle34) > 0.99) && numFreeSurface == numberOfElementNodes && numIsolatedInElement > 1)
                // {
                //     accepted = false;
                //     numberOfSlivers++;
                // }
                // else if ((fabs(cosAngle12) > 0.995 || fabs(cosAngle13) > 0.995 || fabs(cosAngle14) > 0.995 || fabs(cosAngle23) > 0.995 ||
                //           fabs(cosAngle24) > 0.995 || fabs(cosAngle34) > 0.995) && numFreeSurface == numberOfElementNodes && numIsolatedInElement == 1)
                // {
                //     accepted = false;
                //     numberOfSlivers++;
                // }
                if ((fabs(cosAngle12) > 0.999 || fabs(cosAngle13) > 0.999 || fabs(cosAngle14) > 0.999 || fabs(cosAngle23) > 0.999 ||
                          fabs(cosAngle24) > 0.999 || fabs(cosAngle34) > 0.999))
                {
                    accepted = false;
                    numberOfSlivers++;
                }
            }
        }
        if (accepted)
        {
            preservedElements_[el] = ++numberOfAcceptedElements;
        }
    }
    numberOfPreservedElements_ = ++numberOfAcceptedElements;
    info_.numberOfElements_ = numberOfPreservedElements_;

    inMesh_.createPointList(nodes.size(), dimension);

    double* pointList = inMesh_.getPointList();
    int& npoints = inMesh_.getNumberOfPoints();

    int base = 0;
    for (Node* const& node : nodes)
    {
        for (int j = 0; j < dimension; j++)
        {
            double coord = node->getDegreeOfFreedom(j)->getCurrentValue();
            pointList[base+j] = coord;
        }
        base += dimension;
    }
}

void Mesher::generateNewElements(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    int dimension = param->getDimension();

    unsigned int numberOfPrevElements = elements.size(); //number of elements from previous time step
    int* outElementList = outMesh_.getElementList();
    int outNumberOfElements = outMesh_.getNumberOfElements(); //number of elements generated by the mesher
    unsigned int numberOfElementNodes = elements[0]->getNodes().size(); //number of elements that were selected

    //erasing the previous elements
    for (int i = 0; i < numberOfPrevElements; i++)
        delete elements[i]; //after deleting, the pointers in the Geometry class will be invalidated
    elements.clear();
    elements.reserve(numberOfPreservedElements_);

    //reseting isBlocked flag from nodes
    for (Node* const& node : nodes)
        node->setBlocked(false);

    //creating new elements that were selected
    for (int el = 0; el < outNumberOfElements; el++)
    {
        if (preservedElements_[el] >= 0)
        {
            int index = preservedElements_[el];
            std::vector<Node*> elemNodes(numberOfElementNodes);
		    for (int i = 0; i < numberOfElementNodes; i++)
            {
			    elemNodes[i] = nodes[outElementList[el * numberOfElementNodes + i]];
                elemNodes[i]->setBlocked(true);
            }
            Material* mat = nullptr;
            for (Node* const& node : elemNodes)
            {
                if (node->getMaterial())
                {
                    mat = node->getMaterial();  //This should be improved for multi fluid flows
                    break;
                }
            }
            if (!mat) std::cout << "null material detected.\n";

            if (dimension == 2)
            {
                BaseSurfaceElement* base = new BaseSurfaceElement(index, ParametricSurfaceElement::T3, elemNodes);
                base->setPlot(true);
                std::vector<DegreeOfFreedom*> dofs;
                dofs.reserve(elemNodes[0]->getNumberOfDegreesOfFreedom() * numberOfElementNodes);
                for (int i = 0; i < numberOfElementNodes; i++)
                {
                    for (int j = 0; j < dimension; j++)
                    {
                        dofs.emplace_back(elemNodes[i]->getDegreeOfFreedom(j));
                    }
                }
                for (int i = 0; i < numberOfElementNodes; i++)
                {
                    dofs.emplace_back(elemNodes[i]->getDegreeOfFreedom(dimension));
                }
                elements.emplace_back(new PlaneElement(index, dofs, mat, base, param));
            }
            else if (dimension == 3)
            {
                BaseVolumeElement* base = new BaseVolumeElement(index, ParametricVolumeElement::TET4, elemNodes);
                base->setPlot(true);
                std::vector<DegreeOfFreedom*> dofs;
                dofs.reserve(elemNodes[0]->getNumberOfDegreesOfFreedom() * numberOfElementNodes);
                for (int i = 0; i < numberOfElementNodes; i++)
                {
                    for (int j = 0; j < dimension; j++)
                    {
                        dofs.emplace_back(elemNodes[i]->getDegreeOfFreedom(j));
                    }
                }
                for (int i = 0; i < numberOfElementNodes; i++)
                {
                    dofs.emplace_back(elemNodes[i]->getDegreeOfFreedom(dimension));
                }
                elements.emplace_back(new VolumeElement(index, dofs, mat, base, param));
            }
        }
    }

    //Reseting the boundary flags for nodes and elements, except for the nodes that are constrained
    for (Node* const& node : nodes)
        if (!node->isConstrained())
            node->setBoundary(false);
    
    for (Element* const& el : elements)
        el->setBoundary(false);

    //setting the neighbors
    int* outElementNeighborList = outMesh_.getElementNeighbourList();

    int id = 0;
    for (Element* const& el : elements)
    {
        const std::vector<Node*>& elemNodes = el->getNodes();
        for (int i = 0; i < numberOfPreservedElements_; i++)
        {
            if (preservedElements_[id] < 0)
                id++;
            else
                break;
        }

        int numberOfFaces = el->getParametricElement()->getNumberOfFaces();
        int neighborElementId = 0;
        for (int iface = 0; iface < numberOfFaces; iface++)
        {
            neighborElementId = outElementNeighborList[id * numberOfElementNodes + iface];

            if (neighborElementId >= 0 && preservedElements_[neighborElementId] >= 0)
            {
                neighborElementId = preservedElements_[neighborElementId];
                el->addNeighborElement(elements[neighborElementId]);
            }
            else
            {
                el->addNeighborElement(el);
                el->setBoundary(true);
                const std::vector<int>& faceNodes = el->getParametricElement()->getFaceNodes(iface);
                for (int facenode : faceNodes)
                    elemNodes[facenode]->setBoundary(true);
            }
        }
        id++;
    }
}

void Mesher::buildFreeSurface(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    //Reseting the free surface node flag. If isFreeSurface, set the flag of previous time step to true
    for (Node* const& node : nodes)
    {
		if (node->isFreeSurface())
        {
            node->setPreviouslyFreeSurface(true);
            node->setFreeSurface(false);
        }
        else
        {
            node->setPreviouslyFreeSurface(false);
        }
    }

    for (Element* const& el : elements)
    {
        if (el->isBoundary())
        {
            const std::vector<Node*>& elNodes = el->getNodes();
            const std::vector<Element*>& neighborElements = el->getNeighborElements();
            int numberOfFaces = el->getParametricElement()->getNumberOfFaces();

            int face = 0;
            for (Element* const& nElement : neighborElements)
            {
                if (nElement->getIndex() == el->getIndex()) //if the neighbor element is the element itself, the face is a boundary
                {
                    const std::vector<int>& faceVerticesIndex = el->getParametricElement()->getFaceVertices(face);
                    const std::vector<int>& faceNodesIndex = el->getParametricElement()->getFaceNodes(face);
                    unsigned int numberOfVerticesPerFace = faceVerticesIndex.size();
                    unsigned int numberOfNodesPerFace = faceNodesIndex.size();
                    
                    bool freeSurfaceFace = false;
                    for (unsigned int i = 0; i < numberOfVerticesPerFace; i++)
                    {
                        Node* node = elNodes[faceVerticesIndex[i]];
                        if (!node->isConstrained() && !node->isInterface()) //it is sufficient that only one vertice is not an interface nor rigid wall
                        {
                            freeSurfaceFace = true;
                            break;
                        }
                    }
                    if (freeSurfaceFace)
                    {
                        for (unsigned int i = 0; i < numberOfNodesPerFace; i++)
                        {
                            Node* node = elNodes[faceNodesIndex[i]];
                            node->setFreeSurface(true);
                        }
                    }
                }
                face++;
            }
        }
    }
    
}

void Mesher::deactivateSlivers(std::vector<Node*>& nodes, std::vector<Element*>& elements, AnalysisParameters* param)
{
    // This function sets the flag isActive to false for those slivers that are inside the domain.
    // We can't erase them during the remesh procedure, otherwise their nodes would be considered as free surface.
    // Here we simply are not considering their contribution in the system of equations.
    int dimension = param->getDimension();
    double modelVolume = param->getModelVolume();
    double criticalVolume = 0.001 * modelVolume / (double)elements.size();

    if (dimension == 3)
    {
        for (Element* const& el : elements)
	    {
            const std::vector<Node*>& elNodes = el->getNodes();

            int isolatedNodes = 0;
            for (Node* const& node : elNodes)
            {
                const std::vector<Element*>& neighborElements = node->getNeighborElements();
                if (neighborElements.size() == 1)
                    isolatedNodes++;
            }

		    double volume = el->getBaseElement()->getJacobianIntegration();
		    if (volume < criticalVolume && isolatedNodes == 0)
            {
			    el->setActive(false);
            }
            else
            {
                double a1; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                double b1; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                double c1; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,C)
                a1 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) *
                     (elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) - 
                     (elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue());
                b1 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) - 
                     (elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue());
                c1 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue());
                double a2; //slope x for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                double b2; //slope y for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                double c2; //slope z for plane on the second triangular face of the tetrahedra (nodes A,B,D)
                a2 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) *
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue());
                b2 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue());
                c2 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue());
                double a3; //slope x for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                double b3; //slope y for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                double c3; //slope z for plane on the third triangular face of the tetrahedra (nodes B,C,D)
                a3 = (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue());
                b3 = (elNodes[1]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue());
                c3 = (elNodes[1]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[1]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue());
                double a4; //slope x for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                double b4; //slope y for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                double c4; //slope z for plane on the fourth triangular face of the tetrahedra (nodes A,C,D)
                a4 = (elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue());
                b4 = (elNodes[0]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(2)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(2)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue());
                c4 = (elNodes[0]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[3]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue()) -
                     (elNodes[3]->getDegreeOfFreedom(0)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(0)->getCurrentValue()) * 
                     (elNodes[0]->getDegreeOfFreedom(1)->getCurrentValue() - elNodes[2]->getDegreeOfFreedom(1)->getCurrentValue());

                double cosAngle12 = (a1 * a2 + b1 * b2 + c1 * c2) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle13 = (a1 * a3 + b1 * b3 + c1 * c3) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));
                double cosAngle14 = (a1 * a4 + b1 * b4 + c1 * c4) / (sqrt(pow(a1, 2) + pow(b1, 2) + pow(c1, 2)) * sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)));
                double cosAngle23 = (a3 * a2 + b3 * b2 + c3 * c2) / (sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle24 = (a4 * a2 + b4 * b2 + c4 * c2) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a2, 2) + pow(b2, 2) + pow(c2, 2)));
                double cosAngle34 = (a4 * a3 + b4 * b3 + c4 * c3) / (sqrt(pow(a4, 2) + pow(b4, 2) + pow(c4, 2)) * sqrt(pow(a3, 2) + pow(b3, 2) + pow(c3, 2)));

                if ((fabs(cosAngle12) > 0.999 || fabs(cosAngle13) > 0.999 || fabs(cosAngle14) > 0.999 || fabs(cosAngle23) > 0.999 || fabs(cosAngle24) > 0.999 || fabs(cosAngle34) > 0.999) && isolatedNodes == 0)
                {
                    el->setActive(false);
                }
            }
		} 
	}

}
