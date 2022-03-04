#include "VolumeElement.h"
#include <lapacke.h>
#include <algorithm>
#include <iterator>

VolumeElement::VolumeElement(const int& index,
                           const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                           Material* material,
                           BaseVolumeElement* base,
                           AnalysisParameters* params)
        : Element(index, degreesOfFreedons),
          material_(material),
          base_(base),
          parameters_(params) {}

VolumeElement::~VolumeElement()
{
    delete base_;
}

void VolumeElement::setMaterial(Material* material)
{
    material_ = material;
}

void VolumeElement::setBaseElement(BaseVolumeElement* base)
{
    base_ = base;
}

void VolumeElement::setAnalysisParameters(AnalysisParameters* params)
{
    parameters_ = params;
}
        
Material* VolumeElement::getMaterial() const
{
    return material_;
}

BaseVolumeElement* VolumeElement::getBaseElement() const
{
    return base_;
}

ParametricElement* VolumeElement::getParametricElement() const
{
    return base_->getParametricElement();
}

const std::vector<Node*>& VolumeElement::getNodes() const
{
    return base_->getNodes();
}

void VolumeElement::getDOFIndexes(unsigned int& ndof,
                                  int*& indexes) const
{
    ndof = degreesOfFreedom_.size();
    indexes = new int[ndof];

    for (int i = 0; i < ndof; i++)
    {   
        indexes[i] = degreesOfFreedom_[i]->getIndex();
    }
}

void VolumeElement::elementContributions(int& ndofs1,
                                        int& ndofs2,
                                        int*& indexes,
                                        double*& rhsValues,
                                        double*& hessianValues) const
{
    unsigned int ndofs;
    getDOFIndexes(ndofs, indexes);
    unsigned int nterms = ndofs*ndofs;

    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = nodes.size();
    ndofs1 = 3 * nodes.size();              // number of position degrees of freedom
    ndofs2 = ndofs - ndofs1;                // number of pressure degrees of freedom
    bool mixedFormulation = ndofs2;
    
    rhsValues = new double[ndofs];
    hessianValues = new double[nterms];
    
    for (int i = 0; i < nterms; i++)
    {
        hessianValues[i] = 0.0;
    }

    for (int i = 0; i < ndofs; i++)
        rhsValues[i] = 0.0;
    
    const std::vector<QuadraturePoint*>& quadraturePoints = base_->getParametricElement()->getQuadraturePoints();
    const unsigned int numberOfQuadraturePoints = quadraturePoints.size();
    
    double density = material_->getDensity();
    double* gravity = parameters_->getGravity();
    double deltat = parameters_->getDeltat();
    double alphaF = parameters_->getAlphaF();
    double alphaM = parameters_->getAlphaM();
    double gamma = parameters_->getGamma();
    double beta = parameters_->getBeta();

    double tpspg = (0.5 * deltat*deltat) / density; //tpspg = 0.0;

    for (auto& qp : quadraturePoints)
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi[3][3];
        getReferenceJacobianMatrix(dphi_dxsi, dx_dxsi);
        double j0 = getMatrixDeterminant(dx_dxsi);

        double dxsi_dx[3][3];
        getInverseMatrix(dx_dxsi, j0, dxsi_dx);

        double dphi_dx[numberOfNodes][3];
        getReferenceShapeFunctionsGradient(dphi_dxsi, dxsi_dx, dphi_dx);

        double dy_dx[3][3];
        getDeformationGradient(dphi_dx, dy_dx);

        double jacobian = getMatrixDeterminant(dy_dx);

        double dx_dy[3][3];
        getInverseMatrix(dy_dx, jacobian, dx_dy);

        double dphi_dy[numberOfNodes][3];
        for (int i = 0; i < numberOfNodes; i++)
        {
            dphi_dy[i][0] = dphi_dx[i][0] * dx_dy[0][0] + dphi_dx[i][1] * dx_dy[1][0] + dphi_dx[i][2] * dx_dy[2][0];
            dphi_dy[i][1] = dphi_dx[i][0] * dx_dy[0][1] + dphi_dx[i][1] * dx_dy[1][1] + dphi_dx[i][2] * dx_dy[2][1];
            dphi_dy[i][2] = dphi_dx[i][0] * dx_dy[0][2] + dphi_dx[i][1] * dx_dy[1][2] + dphi_dx[i][2] * dx_dy[2][2];
        }

        double CI[6];
        getInverseRightCauchyTensor(dx_dy, CI);
        
        double dE_dy[ndofs1][6];
        getStrainTensorFirstDerivative(dy_dx, dphi_dx, dE_dy);

        double E[6];
        double dv_dx[3][3];
        double dE_dt[6];
        double S[6];
        double dE_dtdy[ndofs1][6];
        double dS_dy[ndofs1][6];

        switch (material_->getType())
        {
            case MaterialType::ELASTIC_SOLID:
                getStrainTensor(dy_dx, E);
                material_->getStressTensorAndDerivative(ndofs1, E, dE_dy, nullptr, S, dS_dy);
                break;
            case MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID:
                getStrainTensor(dy_dx, E);
                material_->getStressTensorAndDerivative(ndofs1, E, dE_dy, nullptr, S, dS_dy);
                break;
            case MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID:
                getDeformationGradientTimeDerivative(dphi_dx, dv_dx);
                getStrainTensorTimeDerivative(dy_dx, dv_dx, dE_dt);
                getStrainTensorTimeDerivativeFirstDerivative(dy_dx, dv_dx, dphi_dx, dE_dtdy);
                material_->getStressTensorAndDerivative(ndofs1, dE_dt, dE_dtdy, dx_dy, S, dS_dy);
                break;
        }

        double accel[3] = {0.0, 0.0, 0.0}; //acceleration at integration point
        for(unsigned int i = 0; i < numberOfNodes; i++)
        {
            accel[0] += phi[i] * nodes[i]->getDegreeOfFreedom(0)->getIntermediateSecondTimeDerivative();
            accel[1] += phi[i] * nodes[i]->getDegreeOfFreedom(1)->getIntermediateSecondTimeDerivative();
            accel[2] += phi[i] * nodes[i]->getDegreeOfFreedom(2)->getIntermediateSecondTimeDerivative();
        }

        double pressure = 0.0; //pressure at integration point
        double dp_dx[3] = {0.0, 0.0, 0.0}; //pressure gradient at integration point
        for(unsigned int i = 0; i < ndofs2; i++)
        {
            pressure += phi[i] * nodes[i]->getDegreeOfFreedom(3)->getCurrentValue();
            dp_dx[0] += dphi_dy[i][0] * nodes[i]->getDegreeOfFreedom(3)->getCurrentValue();
            dp_dx[1] += dphi_dy[i][1] * nodes[i]->getDegreeOfFreedom(3)->getCurrentValue();
            dp_dx[2] += dphi_dy[i][2] * nodes[i]->getDegreeOfFreedom(3)->getCurrentValue();
        }

        const double factor1 = weight * j0;
        const double factor2 = alphaF * factor1;
        const double factor3 = pressure * jacobian;
        const double factor4 = alphaM * density * factor1 / (beta * deltat * deltat);
        const double factor5 = jacobian * factor1;
        const double factor6 = jacobian - 1.0;
        const double factor7 = tpspg * density;
        const double factor8 = tpspg * jacobian;
        const double factor9 = jacobian * factor2;
        const double factor10 = factor8 * factor1;

        {   // Position degrees of freedom
            for (int i = 0; i < ndofs1; i++)
            {
                // node a
                int a = i / 3;
                // dof k
                int k = i % 3;

                // internal force
                double v = doubleContraction(S, dE_dy[i]);

                //inertial and domain force
                double m, bf;
                if (parameters_->useLumpedMass())
                {
                    m = density * nodes[a]->getDegreeOfFreedom(k)->getIntermediateSecondTimeDerivative() / 4.0;
                    bf = density * gravity[k] / 4.0;
                }
                else
                {
                    m = density * phi[a] * accel[k];
                    bf = density * phi[a] * gravity[k];
                }

                //pressure force
                double p = factor3 * doubleContraction(CI, dE_dy[i]);
            
                rhsValues[i] -= (v + m + p - bf) * factor1;
                
                // Position degrees of freedom
                for (int j = 0; j < ndofs1; j++)
                {
                    // node b
                    int b = j / 3;
                    // dof l
                    int l = j % 3;

                    hessianValues[i * ndofs + j] += (doubleContraction(dE_dy[j], dS_dy[i])) * factor2;

                    if (parameters_->useLumpedMass() && i == j)
                        hessianValues[i * ndofs + j] += factor4 / 4.0;

                    if (k == l)
                    {
                        double d2E_dydy[6];
                        getStrainTensorSecondDerivative(a, b, dphi_dx, d2E_dydy);
                        hessianValues[i * ndofs + j] += doubleContraction(S, d2E_dydy) * factor2;
                        if (!parameters_->useLumpedMass())
                            hessianValues[i * ndofs + j] += factor4 * phi[a] * phi[b];
                    }
                }
                // Pressure degrees of freedom
                for (int j = 0; j < ndofs2; j++)
                {
                    // node b
                    int b = j;

                    hessianValues[ndofs1 + i*ndofs + j] += phi[b] * doubleContraction(CI, dE_dy[i]) * factor5;
                }
            }
        }
        {   // Pressure degrees of freedom
            for (int i = 0; i < ndofs2; i++)
            {
                // node a
                int a = i;

                //incompressibility constrain
                double c = phi[a] * factor6;

                //pspg body force
                double bf_pspg = factor7 * (dphi_dy[a][0] * gravity[0] + dphi_dy[a][1] * gravity[1] + dphi_dy[a][2] * gravity[2]); //bf_pspg = 0.0;

                //pspg mass part
                double m_pspg = factor7 * (dphi_dy[a][0] * accel[0] + dphi_dy[a][1] * accel[1] + dphi_dy[a][2] * accel[2]); m_pspg = 0.0;

                //pspg pressure part
                double p_pspg = factor8 * (dphi_dy[a][0] * dp_dx[0] + dphi_dy[a][1] * dp_dx[1] + dphi_dy[a][2] * dp_dx[2]);

                rhsValues[ndofs1 + i] -= (c + m_pspg - p_pspg - bf_pspg) * factor1;
                
                // Position degrees of freedom 
                for (int j = 0; j < ndofs1; j++)
                {
                    // node b
                    int b = j / 3;
                    // dof l
                    int l = j % 3;

                    hessianValues[ndofs*(ndofs1+i) + j] += phi[a] * doubleContraction(CI, dE_dy[j]) * factor9;
                    //hessianValues[ndofs*(ndofs1+i) + j] += tpspg * factor4 * dphi_dy[a][l] * phi[b];
                }
                // Pressure degrees of freedom
                for (int j = 0; j < ndofs2; j++)
                {
                    // node b
                    int b = j;

                    hessianValues[ndofs*(ndofs1+i) + ndofs1 + j] += -(dphi_dy[a][0] * dphi_dy[b][0] + dphi_dy[a][1] * dphi_dy[b][1] + 
                                                                               dphi_dy[a][2] * dphi_dy[b][2]) * factor10;
                }
            }
        }
    }
}

void VolumeElement::clearNeighborElements()
{
    neighborElements_.clear();
    neighborElements_.shrink_to_fit();
    neighborElements_.reserve(base_->getParametricElement()->getNumberOfFaces());
    isBoundary_ = false;
}

void VolumeElement::getCauchyStress(double**& nodalCauchyStress) const
{
    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = nodes.size();
    
    const std::vector<QuadraturePoint*>& quadraturePoints = base_->getParametricElement()->getQuadraturePoints();
    const unsigned int numberOfQuadraturePoints = quadraturePoints.size();
    
    double** gaussCauchyStress = new double*[numberOfQuadraturePoints];
    for (unsigned int i = 0; i < numberOfQuadraturePoints; i++)
        gaussCauchyStress[i] = new double[6];
    
    int sum = -1;
    for (auto& qp : quadraturePoints)
    {
        sum++;
        
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();

        double dx_dxsi[3][3];
        getReferenceJacobianMatrix(dphi_dxsi, dx_dxsi);

        double j0 = getMatrixDeterminant(dx_dxsi);

        double dxsi_dx[3][3];
        getInverseMatrix(dx_dxsi, j0, dxsi_dx);

        double dphi_dx[numberOfNodes][3];
        getReferenceShapeFunctionsGradient(dphi_dxsi, dxsi_dx, dphi_dx);

        double dy_dx[3][3];
        getDeformationGradient(dphi_dx, dy_dx);

        double jacobian = getMatrixDeterminant(dy_dx);

        double dx_dy[3][3];
        getInverseMatrix(dy_dx, jacobian, dx_dy);
        
        double E[6];
        double dv_dx[3][3];
        double dE_dt[6];
        double S[6];

        switch (material_->getType())
        {
            case MaterialType::ELASTIC_SOLID:
                getStrainTensor(dy_dx, E);
                material_->getStressTensor(E, nullptr, S);
                break;
            case MaterialType::ELASTIC_INCOMPRESSIBLE_SOLID:
                getStrainTensor(dy_dx, E);
                material_->getStressTensor(E, nullptr, S);
                break;
            case MaterialType::NEWTONIAN_INCOMPRESSIBLE_FLUID:
                getDeformationGradientTimeDerivative(dphi_dx, dv_dx);
                getStrainTensorTimeDerivative(dy_dx, dv_dx, dE_dt);
                material_->getStressTensor(dE_dt, dx_dy, S);
                break;
        }

        const double factor = 1.0 / jacobian;

        gaussCauchyStress[sum][0] = (dy_dx[0][0] * (dy_dx[0][0] * S[0] + dy_dx[0][1] * S[3] + dy_dx[0][2] * S[4]) +
                                    dy_dx[0][1] * (dy_dx[0][0] * S[3] + dy_dx[0][1] * S[1] + dy_dx[0][2] * S[5]) +
                                    dy_dx[0][2] * (dy_dx[0][0] * S[4] + dy_dx[0][1] * S[5] + dy_dx[0][2] * S[2])) * factor;    // Sigmaxx
        gaussCauchyStress[sum][1] = (dy_dx[1][0] * (dy_dx[1][0] * S[0] + dy_dx[1][1] * S[3] + dy_dx[1][2] * S[4]) +
                                    dy_dx[1][1] * (dy_dx[1][0] * S[3] + dy_dx[1][1] * S[1] + dy_dx[1][2] * S[5]) +
                                    dy_dx[1][2] * (dy_dx[1][0] * S[4] + dy_dx[1][1] * S[5] + dy_dx[1][2] * S[2])) * factor;    // Sigmayy
        gaussCauchyStress[sum][2] = (dy_dx[2][0] * (dy_dx[2][0] * S[0] + dy_dx[2][1] * S[3] + dy_dx[2][2] * S[4]) +
                                    dy_dx[2][1] * (dy_dx[2][0] * S[3] + dy_dx[2][1] * S[1] + dy_dx[2][2] * S[5]) +
                                    dy_dx[2][2] * (dy_dx[2][0] * S[4] + dy_dx[2][1] * S[5] + dy_dx[2][2] * S[2])) * factor;    // Sigmazz
        gaussCauchyStress[sum][3] = (dy_dx[1][0] * (dy_dx[0][0] * S[0] + dy_dx[0][1] * S[3] + dy_dx[0][2] * S[4]) +
                                    dy_dx[1][1] * (dy_dx[0][0] * S[3] + dy_dx[0][1] * S[1] + dy_dx[0][2] * S[5]) +
                                    dy_dx[1][2] * (dy_dx[0][0] * S[4] + dy_dx[0][1] * S[5] + dy_dx[0][2] * S[2])) * factor;    // Sigmaxy
        gaussCauchyStress[sum][4] = (dy_dx[2][0] * (dy_dx[0][0] * S[0] + dy_dx[0][1] * S[3] + dy_dx[0][2] * S[4]) +
                                    dy_dx[2][1] * (dy_dx[0][0] * S[3] + dy_dx[0][1] * S[1] + dy_dx[0][2] * S[5]) +
                                    dy_dx[2][2] * (dy_dx[0][0] * S[4] + dy_dx[0][1] * S[5] + dy_dx[0][2] * S[2])) * factor;    // Sigmaxz
        gaussCauchyStress[sum][5] = (dy_dx[2][0] * (dy_dx[1][0] * S[0] + dy_dx[1][1] * S[3] + dy_dx[1][2] * S[4]) +
                                    dy_dx[2][1] * (dy_dx[1][0] * S[3] + dy_dx[1][1] * S[1] + dy_dx[1][2] * S[5]) +
                                    dy_dx[2][2] * (dy_dx[1][0] * S[4] + dy_dx[1][1] * S[5] + dy_dx[1][2] * S[2])) * factor;    // Sigmayz
    }

    // Transformation matrix
    double* M = new double[numberOfNodes * numberOfNodes];
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        for (unsigned int j = 0; j < numberOfNodes; j++)
        {
            M[numberOfNodes * i + j] = 0.0;
            for (auto& qp : quadraturePoints)
            {
                double* phi = qp->getShapeFunctionsValues();
                M[numberOfNodes * i + j] += phi[i] * phi[j];
            }
        }
    }

    double* cauchyStress = new double[6 * numberOfNodes];
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        cauchyStress[numberOfNodes * 0 + i] = 0.0;
        cauchyStress[numberOfNodes * 1 + i] = 0.0;
        cauchyStress[numberOfNodes * 2 + i] = 0.0;
        cauchyStress[numberOfNodes * 3 + i] = 0.0;
        cauchyStress[numberOfNodes * 4 + i] = 0.0;
        cauchyStress[numberOfNodes * 5 + i] = 0.0;
        for (unsigned int ip = 0; ip < numberOfQuadraturePoints; ip++)
        {
            double* phi = quadraturePoints[ip]->getShapeFunctionsValues();
            cauchyStress[numberOfNodes * 0 + i] += phi[i] * gaussCauchyStress[ip][0];
            cauchyStress[numberOfNodes * 1 + i] += phi[i] * gaussCauchyStress[ip][1];
            cauchyStress[numberOfNodes * 2 + i] += phi[i] * gaussCauchyStress[ip][2];
            cauchyStress[numberOfNodes * 3 + i] += phi[i] * gaussCauchyStress[ip][3];
            cauchyStress[numberOfNodes * 4 + i] += phi[i] * gaussCauchyStress[ip][4];
            cauchyStress[numberOfNodes * 5 + i] += phi[i] * gaussCauchyStress[ip][5];
        }
    }

    int* ipiv = new int[numberOfNodes];
    LAPACKE_dsysv(LAPACK_COL_MAJOR, 'L', numberOfNodes, 6, M, numberOfNodes, ipiv, cauchyStress, numberOfNodes);
    delete[] M;
    delete[] ipiv;

    nodalCauchyStress = new double*[numberOfNodes];
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        nodalCauchyStress[i] = new double[6];
        for (unsigned int j = 0; j < 6; j++)
            nodalCauchyStress[i][j] = cauchyStress[numberOfNodes*j+i];
    }
    if (numberOfQuadraturePoints == 1)
        for (unsigned int i = 0; i < numberOfNodes; i++)
            for (unsigned int j = 0; j < 6; j++)
                nodalCauchyStress[i][j] *= numberOfNodes;

    for (unsigned int ip = 0; ip < numberOfQuadraturePoints; ip++)
    {
        delete[] gaussCauchyStress[ip];
    }
    
    delete[] gaussCauchyStress;
    delete[] cauchyStress;
}

inline void VolumeElement::getReferenceJacobianMatrix(double** dphi_dxsi,
                                                      double A0[3][3]) const
{
    A0[0][0] = 0.0; A0[0][1] = 0.0; A0[0][2] = 0.0;
    A0[1][0] = 0.0; A0[1][1] = 0.0; A0[1][2] = 0.0;
    A0[2][0] = 0.0; A0[2][1] = 0.0; A0[2][2] = 0.0;
    
    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = nodes.size();
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        double x[3];
        switch(referenceConfiguration_)
        {
            case ReferenceConfiguration::INITIAL:
                x[0] = nodes[i]->getDegreeOfFreedom(0)->getInitialValue();
                x[1] = nodes[i]->getDegreeOfFreedom(1)->getInitialValue();
                x[2] = nodes[i]->getDegreeOfFreedom(2)->getInitialValue();
                break;
            case ReferenceConfiguration::PAST:
                x[0] = nodes[i]->getDegreeOfFreedom(0)->getPastValue();
                x[1] = nodes[i]->getDegreeOfFreedom(1)->getPastValue();
                x[2] = nodes[i]->getDegreeOfFreedom(2)->getPastValue();
                break;
            case ReferenceConfiguration::CURRENT:
                std::cout << "Updated Lagrangian Formulation is not implemented on this version of the code.\n";
                exit(EXIT_FAILURE);
        }
        A0[0][0] += dphi_dxsi[0][i] * x[0];
        A0[0][1] += dphi_dxsi[1][i] * x[0];
        A0[0][2] += dphi_dxsi[2][i] * x[0];
        A0[1][0] += dphi_dxsi[0][i] * x[1];
        A0[1][1] += dphi_dxsi[1][i] * x[1];
        A0[1][2] += dphi_dxsi[2][i] * x[1];
        A0[2][0] += dphi_dxsi[0][i] * x[2];
        A0[2][1] += dphi_dxsi[1][i] * x[2];
        A0[2][2] += dphi_dxsi[2][i] * x[2];
    }
}

inline void VolumeElement::getCurrentJacobianMatrix(double** dphi_dxsi,
                                                    double A1[3][3]) const
{
    A1[0][0] = 0.0; A1[0][1] = 0.0; A1[0][2] = 0.0;
    A1[1][0] = 0.0; A1[1][1] = 0.0; A1[1][2] = 0.0;
    A1[2][0] = 0.0; A1[2][1] = 0.0; A1[2][2] = 0.0;
    
    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = nodes.size();
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        double y[3];
        y[0] = nodes[i]->getDegreeOfFreedom(0)->getIntermediateValue();
        y[1] = nodes[i]->getDegreeOfFreedom(1)->getIntermediateValue();
        y[2] = nodes[i]->getDegreeOfFreedom(2)->getIntermediateValue();
        A1[0][0] += dphi_dxsi[0][i] * y[0];
        A1[0][1] += dphi_dxsi[1][i] * y[0];
        A1[0][2] += dphi_dxsi[2][i] * y[0];
        A1[1][0] += dphi_dxsi[0][i] * y[1];
        A1[1][1] += dphi_dxsi[1][i] * y[1];
        A1[1][2] += dphi_dxsi[2][i] * y[1];
        A1[2][0] += dphi_dxsi[0][i] * y[2];
        A1[2][1] += dphi_dxsi[1][i] * y[2];
        A1[2][2] += dphi_dxsi[2][i] * y[2];
    }
}

inline void VolumeElement::getCurrentJacobianMatrixTimeDerivative(double** dphi_dxsi,
                                                                  double dA1_dt[3][3]) const
{
    dA1_dt[0][0] = 0.0; dA1_dt[0][1] = 0.0; dA1_dt[0][2] = 0.0;
    dA1_dt[1][0] = 0.0; dA1_dt[1][1] = 0.0; dA1_dt[1][2] = 0.0;
    dA1_dt[2][0] = 0.0; dA1_dt[2][1] = 0.0; dA1_dt[2][2] = 0.0;
    
    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = nodes.size();
    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        double v[3];
        v[0] = nodes[i]->getDegreeOfFreedom(0)->getIntermediateFirstTimeDerivative();
        v[1] = nodes[i]->getDegreeOfFreedom(1)->getIntermediateFirstTimeDerivative();
        v[2] = nodes[i]->getDegreeOfFreedom(2)->getIntermediateFirstTimeDerivative();
        
        dA1_dt[0][0] += dphi_dxsi[0][i] * v[0];
        dA1_dt[0][1] += dphi_dxsi[1][i] * v[0];
        dA1_dt[0][2] += dphi_dxsi[2][i] * v[0];
        dA1_dt[1][0] += dphi_dxsi[0][i] * v[1];
        dA1_dt[1][1] += dphi_dxsi[1][i] * v[1];
        dA1_dt[1][2] += dphi_dxsi[2][i] * v[1];
        dA1_dt[2][0] += dphi_dxsi[0][i] * v[2];
        dA1_dt[2][1] += dphi_dxsi[1][i] * v[2];
        dA1_dt[2][2] += dphi_dxsi[2][i] * v[2];
    }
}

inline double VolumeElement::getMatrixDeterminant(const double jacobianMatrix[3][3]) const
{
    return jacobianMatrix[0][0] * jacobianMatrix[1][1] * jacobianMatrix[2][2] +
           jacobianMatrix[0][1] * jacobianMatrix[1][2] * jacobianMatrix[2][0] +
           jacobianMatrix[0][2] * jacobianMatrix[1][0] * jacobianMatrix[2][1] -
           jacobianMatrix[0][2] * jacobianMatrix[1][1] * jacobianMatrix[2][0] -
           jacobianMatrix[0][1] * jacobianMatrix[1][0] * jacobianMatrix[2][2] -
           jacobianMatrix[0][0] * jacobianMatrix[1][2] * jacobianMatrix[2][1];
}

inline void VolumeElement::getInverseMatrix(const double matrix[3][3],
                                            const double& determinant,
                                            double inverse[3][3]) const
{
    double inv_determinant = 1.0 / determinant;
    inverse[0][0] = (matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]) * inv_determinant;
    inverse[0][1] = (matrix[0][2] * matrix[2][1] - matrix[0][1] * matrix[2][2]) * inv_determinant;
    inverse[0][2] = (matrix[0][1] * matrix[1][2] - matrix[0][2] * matrix[1][1]) * inv_determinant;
    inverse[1][0] = (matrix[1][2] * matrix[2][0] - matrix[1][0] * matrix[2][2]) * inv_determinant;
    inverse[1][1] = (matrix[0][0] * matrix[2][2] - matrix[0][2] * matrix[2][0]) * inv_determinant;
    inverse[1][2] = (matrix[0][2] * matrix[1][0] - matrix[0][0] * matrix[1][2]) * inv_determinant;
    inverse[2][0] = (matrix[1][0] * matrix[2][1] - matrix[1][1] * matrix[2][0]) * inv_determinant;
    inverse[2][1] = (matrix[0][1] * matrix[2][0] - matrix[0][0] * matrix[2][1]) * inv_determinant;
    inverse[2][2] = (matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]) * inv_determinant;
}

inline void VolumeElement::getReferenceShapeFunctionsGradient(double** dphi_dxsi,
                                                              double dxsi_dx[3][3],
                                                              double dphi_dx[][3]) const
{
    const unsigned int numberOfNodes = base_->getNumberOfNodes();

    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        dphi_dx[i][0] = dphi_dxsi[0][i] * dxsi_dx[0][0] + dphi_dxsi[1][i] * dxsi_dx[1][0] + dphi_dxsi[2][i] * dxsi_dx[2][0];
        dphi_dx[i][1] = dphi_dxsi[0][i] * dxsi_dx[0][1] + dphi_dxsi[1][i] * dxsi_dx[1][1] + dphi_dxsi[2][i] * dxsi_dx[2][1];
        dphi_dx[i][2] = dphi_dxsi[0][i] * dxsi_dx[0][2] + dphi_dxsi[1][i] * dxsi_dx[1][2] + dphi_dxsi[2][i] * dxsi_dx[2][2];
    }
}

inline void VolumeElement::getCurrentShapeFunctionsGradient(double** dphi_dxsi,
                                                            double dxsi_dy[3][3],
                                                            double dphi_dy[][3]) const
{
    const unsigned int numberOfNodes = base_->getNumberOfNodes();

    for (unsigned int i = 0; i < numberOfNodes; i++)
    {
        dphi_dy[i][0] = dphi_dxsi[0][i] * dxsi_dy[0][0] + dphi_dxsi[1][i] * dxsi_dy[1][0] + dphi_dxsi[2][i] * dxsi_dy[2][0];
        dphi_dy[i][1] = dphi_dxsi[0][i] * dxsi_dy[0][1] + dphi_dxsi[1][i] * dxsi_dy[1][1] + dphi_dxsi[2][i] * dxsi_dy[2][1];
        dphi_dy[i][2] = dphi_dxsi[0][i] * dxsi_dy[0][2] + dphi_dxsi[1][i] * dxsi_dy[1][2] + dphi_dxsi[2][i] * dxsi_dy[2][2];
    }
}

inline void VolumeElement::getDeformationGradient(double dphi_dx[][3],
                                                  double A[3][3]) const
{
    A[0][0] = 0.0; A[0][1] = 0.0; A[0][2] = 0.0;
    A[1][0] = 0.0; A[1][1] = 0.0; A[1][2] = 0.0;
    A[2][0] = 0.0; A[2][1] = 0.0; A[2][2] = 0.0;

    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = base_->getNumberOfNodes();
    for (int i = 0; i < numberOfNodes; i++)
    {
        double y[3];
        y[0] = nodes[i]->getDegreeOfFreedom(0)->getIntermediateValue();
        y[1] = nodes[i]->getDegreeOfFreedom(1)->getIntermediateValue();
        y[2] = nodes[i]->getDegreeOfFreedom(2)->getIntermediateValue();
        A[0][0] += dphi_dx[i][0] * y[0];
        A[0][1] += dphi_dx[i][1] * y[0];
        A[0][2] += dphi_dx[i][2] * y[0];
        A[1][0] += dphi_dx[i][0] * y[1];
        A[1][1] += dphi_dx[i][1] * y[1];
        A[1][2] += dphi_dx[i][2] * y[1];
        A[2][0] += dphi_dx[i][0] * y[2];
        A[2][1] += dphi_dx[i][1] * y[2];
        A[2][2] += dphi_dx[i][2] * y[2];
    }
}

inline void VolumeElement::getDeformationGradientTimeDerivative(double dphi_dx[][3],
                                                                double dA_dt[3][3]) const
{
    dA_dt[0][0] = 0.0; dA_dt[0][1] = 0.0; dA_dt[0][2] = 0.0;
    dA_dt[1][0] = 0.0; dA_dt[1][1] = 0.0; dA_dt[1][2] = 0.0;
    dA_dt[2][0] = 0.0; dA_dt[2][1] = 0.0; dA_dt[2][2] = 0.0;

    const std::vector<Node*>& nodes = base_->getNodes();
    const unsigned int numberOfNodes = base_->getNumberOfNodes();

    for (int i = 0; i < numberOfNodes; i++)
    {
        double v[3];
        v[0] = nodes[i]->getDegreeOfFreedom(0)->getIntermediateFirstTimeDerivative();
        v[1] = nodes[i]->getDegreeOfFreedom(1)->getIntermediateFirstTimeDerivative();
        v[2] = nodes[i]->getDegreeOfFreedom(2)->getIntermediateFirstTimeDerivative();
        dA_dt[0][0] += dphi_dx[i][0] * v[0];
        dA_dt[0][1] += dphi_dx[i][1] * v[0];
        dA_dt[0][2] += dphi_dx[i][2] * v[0];
        dA_dt[1][0] += dphi_dx[i][0] * v[1];
        dA_dt[1][1] += dphi_dx[i][1] * v[1];
        dA_dt[1][2] += dphi_dx[i][2] * v[1];
        dA_dt[2][0] += dphi_dx[i][0] * v[2];
        dA_dt[2][1] += dphi_dx[i][1] * v[2];
        dA_dt[2][2] += dphi_dx[i][2] * v[2];
    }
}

inline void VolumeElement::getStrainTensor(const double A[3][3],
                                           double E[6]) const
{
    E[0] = 0.5 * (A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0]) - 0.5;     //E[0][0]
    E[1] = 0.5 * (A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]) - 0.5;     //E[1][1]
    E[2] = 0.5 * (A[0][2] * A[0][2] + A[1][2] * A[1][2] + A[2][2] * A[2][2]) - 0.5;     //E[2][2]
    E[3] = 0.5 * (A[0][0] * A[0][1] + A[1][0] * A[1][1] + A[2][0] * A[2][1]);           //E[0][1]
    E[4] = 0.5 * (A[0][0] * A[0][2] + A[1][0] * A[1][2] + A[2][0] * A[2][2]);           //E[0][2]
    E[5] = 0.5 * (A[0][1] * A[0][2] + A[1][1] * A[1][2] + A[2][1] * A[2][2]);           //E[1][2]
}

inline void VolumeElement::getStrainTensorTimeDerivative(const double A[3][3],
                                                         const double dA_dt[3][3],
                                                         double dE_dt[6]) const
{
    dE_dt[0] = A[0][0] * dA_dt[0][0] + A[1][0] * dA_dt[1][0] + A[2][0] * dA_dt[2][0];   //dE_dt[0][0]
    dE_dt[1] = A[0][1] * dA_dt[0][1] + A[1][1] * dA_dt[1][1] + A[2][1] * dA_dt[2][1];   //dE_dt[1][1]
    dE_dt[2] = A[0][2] * dA_dt[0][2] + A[1][2] * dA_dt[1][2] + A[2][2] * dA_dt[2][2];   //dE_dt[2][2]
    dE_dt[3] = 0.5 * (dA_dt[0][0] * A[0][1] + A[0][0] * dA_dt[0][1] +
                      dA_dt[1][0] * A[1][1] + A[1][0] * dA_dt[1][1] +
                      dA_dt[2][0] * A[2][1] + A[2][0] * dA_dt[2][1]);       //dE_dt[0][1]
    dE_dt[4] = 0.5 * (dA_dt[0][0] * A[0][2] + A[0][0] * dA_dt[0][2] +
                      dA_dt[1][0] * A[1][2] + A[1][0] * dA_dt[1][2] +
                      dA_dt[2][0] * A[2][2] + A[2][0] * dA_dt[2][2]);       //dE_dt[0][2]
    dE_dt[5] = 0.5 * (dA_dt[0][1] * A[0][2] + A[0][1] * dA_dt[0][2] +
                      dA_dt[1][1] * A[1][2] + A[1][1] * dA_dt[1][2] +
                      dA_dt[2][1] * A[2][2] + A[2][1] * dA_dt[2][2]);       //dE_dt[1][2]
}

inline void VolumeElement::getInverseRightCauchyTensor(const double AI[3][3],
                                                       double CI[6]) const
{
    CI[0] = AI[0][0] * AI[0][0] + AI[0][1] * AI[0][1] + AI[0][2] * AI[0][2];    //CI[0][0]
    CI[1] = AI[1][0] * AI[1][0] + AI[1][1] * AI[1][1] + AI[1][2] * AI[1][2];    //CI[1][1]
    CI[2] = AI[2][0] * AI[2][0] + AI[2][1] * AI[2][1] + AI[2][2] * AI[2][2];    //CI[2][2]
    CI[3] = AI[0][0] * AI[1][0] + AI[0][1] * AI[1][1] + AI[0][2] * AI[1][2];    //CI[0][1]
    CI[4] = AI[0][0] * AI[2][0] + AI[0][1] * AI[2][1] + AI[0][2] * AI[2][2];    //CI[0][2]
    CI[5] = AI[1][0] * AI[2][0] + AI[1][1] * AI[2][1] + AI[1][2] * AI[2][2];    //CI[1][2]
}

inline void VolumeElement::getStrainTensorFirstDerivative(const double A[3][3],
                                                          double dphi_dx[][3],
                                                          double dE_dy[][6]) const
{
    const unsigned int numberOfNodes = base_->getNumberOfNodes();
    const unsigned int nPositionDOFs = 3 * numberOfNodes;

    for (unsigned int i = 0; i < nPositionDOFs; i++)
    {
        const int a = i / 3; //node a
        const int j = i % 3; //dof j

        dE_dy[i][0] = dphi_dx[a][0] * A[j][0];                                    //dE_dy[0][0]
        dE_dy[i][1] = dphi_dx[a][1] * A[j][1];                                    //dE_dy[1][1]
        dE_dy[i][2] = dphi_dx[a][2] * A[j][2];                                    //dE_dy[2][2]
        dE_dy[i][3] = 0.5 * (dphi_dx[a][0] * A[j][1] + A[j][0] * dphi_dx[a][1]);  //dE_dy[0][1]
        dE_dy[i][4] = 0.5 * (dphi_dx[a][0] * A[j][2] + A[j][0] * dphi_dx[a][2]);  //dE_dy[0][2]
        dE_dy[i][5] = 0.5 * (dphi_dx[a][1] * A[j][2] + A[j][1] * dphi_dx[a][2]);  //dE_dy[1][2]
    }
}

inline void VolumeElement::getStrainTensorTimeDerivativeFirstDerivative(const double A[3][3],
                                                                        const double dA_dt[3][3],
                                                                        double dphi_dx[][3],
                                                                        double dE_dtdy[][6]) const
{
    const unsigned int numberOfNodes = base_->getNumberOfNodes();
    const unsigned int nPositionDOFs = 3 * numberOfNodes;

    const double factor = parameters_->getGamma() / (parameters_->getBeta() * parameters_->getDeltat());

    for (unsigned int i = 0; i < nPositionDOFs; i++)
    {
        int a = i / 3; //node a
        int j = i % 3; //dof j

        dE_dtdy[i][0] = factor * dphi_dx[a][0] * A[j][0] + dA_dt[j][0] * dphi_dx[a][0];             //dE_dtdy[0][0]
        dE_dtdy[i][1] = factor * dphi_dx[a][1] * A[j][1] + dA_dt[j][1] * dphi_dx[a][1];             //dE_dtdy[1][1]
        dE_dtdy[i][2] = factor * dphi_dx[a][2] * A[j][2] + dA_dt[j][2] * dphi_dx[a][2];             //dE_dtdy[2][2]
        dE_dtdy[i][3] = 0.5 * (factor * dphi_dx[a][0] * A[j][1] + dA_dt[j][0] * dphi_dx[a][1] +
                               factor * dphi_dx[a][1] * A[j][0] + dA_dt[j][1] * dphi_dx[a][0]);     //dE_dtdy[0][1]
        dE_dtdy[i][4] = 0.5 * (factor * dphi_dx[a][0] * A[j][2] + dA_dt[j][0] * dphi_dx[a][2] +
                               factor * dphi_dx[a][2] * A[j][0] + dA_dt[j][2] * dphi_dx[a][0]);     //dE_dtdy[0][2]
        dE_dtdy[i][5] = 0.5 * (factor * dphi_dx[a][1] * A[j][2] + dA_dt[j][1] * dphi_dx[a][2] +
                               factor * dphi_dx[a][2] * A[j][1] + dA_dt[j][2] * dphi_dx[a][1]);     //dE_dtdy[1][2]
    }
}

inline void VolumeElement::getStrainTensorSecondDerivative(const int i,
                                                           const int j,
                                                           const double dphi_dx[][3],
                                                           double d2E_dydy[6]) const
{
    d2E_dydy[0] =  dphi_dx[i][0] * dphi_dx[j][0];                                        // d2E_dydy[0][0]
    d2E_dydy[1] =  dphi_dx[i][1] * dphi_dx[j][1];                                        // d2E_dydy[1][1]
    d2E_dydy[2] =  dphi_dx[i][2] * dphi_dx[j][2];                                        // d2E_dydy[2][2]
    d2E_dydy[3] = 0.5 * (dphi_dx[i][0] * dphi_dx[j][1] + dphi_dx[j][0] * dphi_dx[i][1]); // d2E_dydy[0][1]
    d2E_dydy[4] = 0.5 * (dphi_dx[i][0] * dphi_dx[j][2] + dphi_dx[j][0] * dphi_dx[i][2]); // d2E_dydy[0][2]
    d2E_dydy[5] = 0.5 * (dphi_dx[i][1] * dphi_dx[j][2] + dphi_dx[j][1] * dphi_dx[i][2]); // d2E_dydy[1][2]
}

inline double VolumeElement::doubleContraction(const double v1[6],
                                               const double v2[6]) const
{
    return  v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] +
           (v1[3] * v2[3] + v1[4] * v2[4] + v1[5] * v2[5]) * 2.0;
}