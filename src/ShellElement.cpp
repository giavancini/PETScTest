#include "ShellElement.h"
#include "Quadratures.h"
#include <lapacke.h>

ShellElement::ShellElement(const int& index,
                           const std::vector<DegreeOfFreedom*>& degreesOfFreedons,
                           const double& thickness,
                           Material* material,
                           BaseSurfaceElement* base,
                           AnalysisParameters* params)
    : Element(index, degreesOfFreedons),
      thickness_(thickness),
      material_(material),
      base_(base),
      parameters_(params),
      numberOfThicknessIntegrationPoints_(3) {}

ShellElement::~ShellElement()
{
    delete base_;
}

void ShellElement::setMaterial(Material* material)
{
    material_ = material;
}

void ShellElement::setBaseElement(BaseSurfaceElement* base)
{
    base_ = base;
}

void ShellElement::setAnalysisParameters(AnalysisParameters* params)
{
    parameters_ = params;
}
        
Material* ShellElement::getMaterial() const
{
    return material_;
}

BaseSurfaceElement* ShellElement::getBaseElement() const
{
    return base_;
}

const std::vector<Node*>& ShellElement::getNodes() const
{
    return base_->getNodes();
}

void ShellElement::getDOFIndexes(unsigned int& ndof,
                                  int*& indexes) const
{
    ndof = degreesOfFreedom_.size();
    indexes = new int[ndof];

    for (int i = 0; i < ndof; i++)
    {   
        indexes[i] = degreesOfFreedom_[i]->getIndex();
    }
}

void ShellElement::elementContributions(int& ndofs1,
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
    ndofs1 = 7 * nodes.size();              // number of position, generalized vector and thickness deformation rate degrees of freedom
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

    // surface parameters
    const std::vector<QuadraturePoint*>& quadraturePoints = base_->getParametricElement()->getQuadraturePoints();
    const unsigned int numberOfQuadraturePoints = quadraturePoints.size();

    // thickness parameters
    double** xsi3;
    double* weightThickness;
    quadratures::lineQuadrature(numberOfThicknessIntegrationPoints_, xsi3, weightThickness);

    double density = material_->getDensity();
    double* gravity = parameters_->getGravity();
    double deltat = parameters_->getDeltat();
    double alphaF = parameters_->getAlphaF();
    double alphaM = parameters_->getAlphaM();
    double gamma = parameters_->getGamma();
    double beta = parameters_->getBeta();

    for (auto& qp : quadraturePoints) //Surface quadrature points
    {
        double* phi = qp->getShapeFunctionsValues();
        double** dphi_dxsi = qp->getShapeFunctionsDerivativesValues();
        double weight = qp->getWeight();
        
        double n[3], g[3], t0, t, dx_dxsi[3][2], dy_dxsi[3][2], dn_dxsi[3][2], dg_dxsi[3][2], dt0_dxsi[2], dt_dxsi[2];
        getSurfaceDependentVariables(phi, dphi_dxsi, n, g, t0, t, dx_dxsi, dy_dxsi, dn_dxsi, dg_dxsi, dt0_dxsi, dt_dxsi);

        double accel[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //acceleration at integration point
        for(unsigned int i = 0; i < numberOfNodes; i++)
        {
            accel[0] += phi[i] * nodes[i]->getDegreeOfFreedom(0)->getIntermediateSecondTimeDerivative();
            accel[1] += phi[i] * nodes[i]->getDegreeOfFreedom(1)->getIntermediateSecondTimeDerivative();
            accel[2] += phi[i] * nodes[i]->getDegreeOfFreedom(2)->getIntermediateSecondTimeDerivative();
            accel[3] += phi[i] * nodes[i]->getDegreeOfFreedom(3)->getIntermediateSecondTimeDerivative();
            accel[4] += phi[i] * nodes[i]->getDegreeOfFreedom(4)->getIntermediateSecondTimeDerivative();
            accel[5] += phi[i] * nodes[i]->getDegreeOfFreedom(5)->getIntermediateSecondTimeDerivative();
        }

        for (int ig = 0; ig < numberOfThicknessIntegrationPoints_; ig++)
        {
            const double hm = 0.5 * thickness_ * xsi3[0][ig];

            double A0[3][3];
            getReferenceJacobianMatrix(t0, thickness_, hm, n, dx_dxsi, dn_dxsi, dt0_dxsi, A0);

            double j0 = getMatrixDeterminant(A0);

            double A1[3][3];
            getCurrentJacobianMatrix(t, thickness_, hm, g, dy_dxsi, dg_dxsi, dt_dxsi, A1);

            double A0I[3][3];
            getInverseMatrix(A0, j0, A0I);

            double A[3][3];
            getDeformationGradient(A1, A0I, A);

            double E[6];
            getStrainTensor(A, E);

            double dA_dy[ndofs1][3][3];
            double dE_dy[ndofs1][6];
            getFirstDerivatives(phi, dphi_dxsi, hm, thickness_, t, g, dt_dxsi, A0I, A, dg_dxsi, dA_dy, dE_dy);
        
            double S[6];
            double dS_dy[ndofs1][6];
            material_->getStressTensorAndDerivative(ndofs1, E, dE_dy, nullptr, S, dS_dy);

            double accel2[3];
            accel2[0] = accel[0] + hm * accel[3];
            accel2[1] = accel[1] + hm * accel[4];
            accel2[2] = accel[2] + hm * accel[5];

            const double factor1 = weight * weightThickness[ig] * j0;
            const double factor2 = alphaF * factor1;
            const double factor3 = alphaM * density * factor1 / (beta * deltat * deltat);

            for (int i = 0; i < ndofs1; i++)
            {
                // node a
                int a = i / 7;
                // dof k
                int k = i % 7;

                int dof_mass = (k < 3) ? k : k - 3;

                // internal force
                double v = doubleContraction(S, dE_dy[i]);

                //inertial force
                double m = (k != 6) ? density * phi[a] * accel2[dof_mass] : 0.0;

                //domain force
                double bf = density * phi[a] * gravity[k];
            
                rhsValues[i] -= (v + m - bf) * factor1;
                
                // Position degrees of freedom
                for (int j = 0; j < ndofs1; j++)
                {
                    // node b
                    int b = j / 7;
                    // dof l
                    int l = j % 7;

                    double d2E_dydy[6];
                    getStrainTensorSecondDerivative(phi, dphi_dxsi, thickness_, hm, i, j, A0I, A, dA_dy, d2E_dydy);

                    hessianValues[i * ndofs + j] += (doubleContraction(dE_dy[j], dS_dy[i]) + doubleContraction(S, d2E_dydy)) * factor2;

                    if (k == 6 || l == 6)
                    {
                        continue;
                    }
                    else if (k < 3 && k == l)
                    {
                        hessianValues[i * ndofs + j] += factor3 * phi[a] * phi[b];
                    }
                    else if (k == l + 3 || l == k + 3)
                    {
                        hessianValues[i * ndofs + j] += factor3 * hm * phi[a] * phi[b];
                    }
                    else if (k >= 3 && k == l)
                    {
                        hessianValues[i * ndofs + j] += factor3 * hm * hm * phi[a] * phi[b];
                    }
                }
            }
        
        }
    }
}

void ShellElement::getCauchyStress(double**& nodalCauchyStress) const
{

}

void ShellElement::incrementNormalVector()
{
    const int numberOfNodes = base_->getNumberOfNodes();
    const std::vector<Node*>& nodes = getNodes();
    double** nodalParametricCoordinates = base_->getParametricElement()->getNodalParametricCoordinates();
    for (int i = 0; i < numberOfNodes; i++)
    {
        double normal[3];
        double xsi[2] = {nodalParametricCoordinates[0][i], nodalParametricCoordinates[1][i]}; 
        base_->getInitialNormalVector(xsi, normal);

        for (int j = 0; j < 3; j++)
        {
            double value = nodes[i]->getDegreeOfFreedom(3+j)->getInitialValue() + normal[j];
            nodes[i]->getDegreeOfFreedom(3+j)->setInitialValue(value);
        }
    }
}

inline void ShellElement::getSurfaceDependentVariables(double* phi,
                                                       double** dphi_dxsi,
                                                       double n[3],
                                                       double g[3],
                                                       double& t0,
                                                       double& t,
                                                       double dx_dxsi[3][2],
                                                       double dy_dxsi[3][2],
                                                       double dn_dxsi[3][2],
                                                       double dg_dxsi[3][2],
                                                       double dt0_dxsi[2],
                                                       double dt_dxsi[2]) const
{
    n[0] = 0.0; n[1] = 0.0; n[2] = 0.0;
    g[0] = 0.0; g[1] = 0.0; g[2] = 0.0;
    t0 = 0.0;
    t = 0.0;
    dx_dxsi[0][0] = 0.0; dx_dxsi[0][1] = 0.0; dx_dxsi[1][0] = 0.0;
    dx_dxsi[1][1] = 0.0; dx_dxsi[2][0] = 0.0; dx_dxsi[2][1] = 0.0;
    dy_dxsi[0][0] = 0.0; dy_dxsi[0][1] = 0.0; dy_dxsi[1][0] = 0.0;
    dy_dxsi[1][1] = 0.0; dy_dxsi[2][0] = 0.0; dy_dxsi[2][1] = 0.0;
    dn_dxsi[0][0] = 0.0; dn_dxsi[0][1] = 0.0; dn_dxsi[1][0] = 0.0;
    dn_dxsi[1][1] = 0.0; dn_dxsi[2][0] = 0.0; dn_dxsi[2][1] = 0.0;
    dg_dxsi[0][0] = 0.0; dg_dxsi[0][1] = 0.0; dg_dxsi[1][0] = 0.0;
    dg_dxsi[1][1] = 0.0; dg_dxsi[2][0] = 0.0; dg_dxsi[2][1] = 0.0;
    dt0_dxsi[0] = 0.0; dt0_dxsi[1] = 0.0;
    dt_dxsi[0] = 0.0; dt_dxsi[1] = 0.0;

    const int numberOfNodes = base_->getNumberOfNodes();
    const std::vector<Node*>& nodes = getNodes();
    for (int i = 0; i < numberOfNodes; ++i)
    {
        n[0] += nodes[i]->getDegreeOfFreedom(3)->getInitialValue() * phi[i];
        n[1] += nodes[i]->getDegreeOfFreedom(4)->getInitialValue() * phi[i];
        n[2] += nodes[i]->getDegreeOfFreedom(5)->getInitialValue() * phi[i];

        g[0] += nodes[i]->getDegreeOfFreedom(3)->getIntermediateValue() * phi[i];
        g[1] += nodes[i]->getDegreeOfFreedom(4)->getIntermediateValue() * phi[i];
        g[2] += nodes[i]->getDegreeOfFreedom(5)->getIntermediateValue() * phi[i];

        t0 +=  nodes[i]->getDegreeOfFreedom(6)->getInitialValue() * phi[i];

        t +=  nodes[i]->getDegreeOfFreedom(6)->getIntermediateValue() * phi[i];

        dx_dxsi[0][0] += nodes[i]->getDegreeOfFreedom(0)->getInitialValue() * dphi_dxsi[0][i];
        dx_dxsi[0][1] += nodes[i]->getDegreeOfFreedom(0)->getInitialValue() * dphi_dxsi[1][i];
        dx_dxsi[1][0] += nodes[i]->getDegreeOfFreedom(1)->getInitialValue() * dphi_dxsi[0][i];
        dx_dxsi[1][1] += nodes[i]->getDegreeOfFreedom(1)->getInitialValue() * dphi_dxsi[1][i];
        dx_dxsi[2][0] += nodes[i]->getDegreeOfFreedom(2)->getInitialValue() * dphi_dxsi[0][i];
        dx_dxsi[2][1] += nodes[i]->getDegreeOfFreedom(2)->getInitialValue() * dphi_dxsi[1][i];

        dy_dxsi[0][0] += nodes[i]->getDegreeOfFreedom(0)->getIntermediateValue() * dphi_dxsi[0][i];
        dy_dxsi[0][1] += nodes[i]->getDegreeOfFreedom(0)->getIntermediateValue() * dphi_dxsi[1][i];
        dy_dxsi[1][0] += nodes[i]->getDegreeOfFreedom(1)->getIntermediateValue() * dphi_dxsi[0][i];
        dy_dxsi[1][1] += nodes[i]->getDegreeOfFreedom(1)->getIntermediateValue() * dphi_dxsi[1][i];
        dy_dxsi[2][0] += nodes[i]->getDegreeOfFreedom(2)->getIntermediateValue() * dphi_dxsi[0][i];
        dy_dxsi[2][1] += nodes[i]->getDegreeOfFreedom(2)->getIntermediateValue() * dphi_dxsi[1][i];

        dn_dxsi[0][0] += nodes[i]->getDegreeOfFreedom(3)->getInitialValue() * dphi_dxsi[0][i];
        dn_dxsi[0][1] += nodes[i]->getDegreeOfFreedom(3)->getInitialValue() * dphi_dxsi[1][i];
        dn_dxsi[1][0] += nodes[i]->getDegreeOfFreedom(4)->getInitialValue() * dphi_dxsi[0][i];
        dn_dxsi[1][1] += nodes[i]->getDegreeOfFreedom(4)->getInitialValue() * dphi_dxsi[1][i];
        dn_dxsi[2][0] += nodes[i]->getDegreeOfFreedom(5)->getInitialValue() * dphi_dxsi[0][i];
        dn_dxsi[2][1] += nodes[i]->getDegreeOfFreedom(5)->getInitialValue() * dphi_dxsi[1][i];

        dg_dxsi[0][0] += nodes[i]->getDegreeOfFreedom(3)->getIntermediateValue() * dphi_dxsi[0][i];
        dg_dxsi[0][1] += nodes[i]->getDegreeOfFreedom(3)->getIntermediateValue() * dphi_dxsi[1][i];
        dg_dxsi[1][0] += nodes[i]->getDegreeOfFreedom(4)->getIntermediateValue() * dphi_dxsi[0][i];
        dg_dxsi[1][1] += nodes[i]->getDegreeOfFreedom(4)->getIntermediateValue() * dphi_dxsi[1][i];
        dg_dxsi[2][0] += nodes[i]->getDegreeOfFreedom(5)->getIntermediateValue() * dphi_dxsi[0][i];
        dg_dxsi[2][1] += nodes[i]->getDegreeOfFreedom(5)->getIntermediateValue() * dphi_dxsi[1][i];

        dt0_dxsi[0] += nodes[i]->getDegreeOfFreedom(6)->getInitialValue() * dphi_dxsi[0][i];
        dt0_dxsi[1] += nodes[i]->getDegreeOfFreedom(6)->getInitialValue() * dphi_dxsi[1][i];

        dt_dxsi[0] += nodes[i]->getDegreeOfFreedom(6)->getIntermediateValue() * dphi_dxsi[0][i];
        dt_dxsi[1] += nodes[i]->getDegreeOfFreedom(6)->getIntermediateValue() * dphi_dxsi[1][i];
    }
}

inline void ShellElement::getReferenceJacobianMatrix(const double t0,
                                                     const double h0,
                                                     const double hm,
                                                     const double n[3],
                                                     const double dx_dxsi[3][2],
                                                     const double dn_dxsi[3][2],
                                                     const double dt0_dxsi[2],
                                                     double A0[3][3]) const
{
    A0[0][0] = dx_dxsi[0][0] + hm * ((1.0 + t0 * hm) * dn_dxsi[0][0] + dt0_dxsi[0] * n[0] * hm);
    A0[0][1] = dx_dxsi[0][1] + hm * ((1.0 + t0 * hm) * dn_dxsi[0][1] + dt0_dxsi[1] * n[0] * hm);
    A0[0][2] = 0.5 * h0 * (1.0 + 2.0 * t0 * hm) * n[0];
    A0[1][0] = dx_dxsi[1][0] + hm * ((1.0 + t0 * hm) * dn_dxsi[1][0] + dt0_dxsi[0] * n[1] * hm);
    A0[1][1] = dx_dxsi[1][1] + hm * ((1.0 + t0 * hm) * dn_dxsi[1][1] + dt0_dxsi[1] * n[1] * hm);
    A0[1][2] = 0.5 * h0 * (1.0 + 2.0 * t0 * hm) * n[1];
    A0[2][0] = dx_dxsi[2][0] + hm * ((1.0 + t0 * hm) * dn_dxsi[2][0] + dt0_dxsi[0] * n[2] * hm);
    A0[2][1] = dx_dxsi[2][1] + hm * ((1.0 + t0 * hm) * dn_dxsi[2][1] + dt0_dxsi[1] * n[2] * hm);
    A0[2][2] = 0.5 * h0 * (1.0 + 2.0 * t0 * hm) * n[2];
}

inline void ShellElement::getCurrentJacobianMatrix(const double t,
                                                   const double h0,
                                                   const double hm,
                                                   const double g[3],
                                                   const double dy_dxsi[3][2],
                                                   const double dg_dxsi[3][2],
                                                   const double dt_dxsi[2],
                                                   double A1[3][3]) const
{
    A1[0][0] = dy_dxsi[0][0] + hm * ((1.0 + t * hm) * dg_dxsi[0][0] + dt_dxsi[0] * g[0] * hm);
    A1[0][1] = dy_dxsi[0][1] + hm * ((1.0 + t * hm) * dg_dxsi[0][1] + dt_dxsi[1] * g[0] * hm);
    A1[0][2] = 0.5 * h0 * (1.0 + 2.0 * t * hm) * g[0];
    A1[1][0] = dy_dxsi[1][0] + hm * ((1.0 + t * hm) * dg_dxsi[1][0] + dt_dxsi[0] * g[1] * hm);
    A1[1][1] = dy_dxsi[1][1] + hm * ((1.0 + t * hm) * dg_dxsi[1][1] + dt_dxsi[1] * g[1] * hm);
    A1[1][2] = 0.5 * h0 * (1.0 + 2.0 * t * hm) * g[1];
    A1[2][0] = dy_dxsi[2][0] + hm * ((1.0 + t * hm) * dg_dxsi[2][0] + dt_dxsi[0] * g[2] * hm);
    A1[2][1] = dy_dxsi[2][1] + hm * ((1.0 + t * hm) * dg_dxsi[2][1] + dt_dxsi[1] * g[2] * hm);
    A1[2][2] = 0.5 * h0 * (1.0 + 2.0 * t * hm) * g[2];
}

inline double ShellElement::getMatrixDeterminant(const double jacobianMatrix[3][3]) const
{
    return jacobianMatrix[0][0] * jacobianMatrix[1][1] * jacobianMatrix[2][2] +
           jacobianMatrix[0][1] * jacobianMatrix[1][2] * jacobianMatrix[2][0] +
           jacobianMatrix[0][2] * jacobianMatrix[1][0] * jacobianMatrix[2][1] -
           jacobianMatrix[0][2] * jacobianMatrix[1][1] * jacobianMatrix[2][0] -
           jacobianMatrix[0][1] * jacobianMatrix[1][0] * jacobianMatrix[2][2] -
           jacobianMatrix[0][0] * jacobianMatrix[1][2] * jacobianMatrix[2][1];
}

inline void ShellElement::getInverseMatrix(const double matrix[3][3],
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

inline void ShellElement::getDeformationGradient(const double A1[3][3],
                                                 const double A0I[3][3],
                                                 double A[3][3]) const
{
    A[0][0] = A1[0][0] * A0I[0][0] + A1[0][1] * A0I[1][0] + A1[0][2] * A0I[2][0];
    A[0][1] = A1[0][0] * A0I[0][1] + A1[0][1] * A0I[1][1] + A1[0][2] * A0I[2][1];
    A[0][2] = A1[0][0] * A0I[0][2] + A1[0][1] * A0I[1][2] + A1[0][2] * A0I[2][2];
    A[1][0] = A1[1][0] * A0I[0][0] + A1[1][1] * A0I[1][0] + A1[1][2] * A0I[2][0];
    A[1][1] = A1[1][0] * A0I[0][1] + A1[1][1] * A0I[1][1] + A1[1][2] * A0I[2][1];
    A[1][2] = A1[1][0] * A0I[0][2] + A1[1][1] * A0I[1][2] + A1[1][2] * A0I[2][2];
    A[2][0] = A1[2][0] * A0I[0][0] + A1[2][1] * A0I[1][0] + A1[2][2] * A0I[2][0];
    A[2][1] = A1[2][0] * A0I[0][1] + A1[2][1] * A0I[1][1] + A1[2][2] * A0I[2][1];
    A[2][2] = A1[2][0] * A0I[0][2] + A1[2][1] * A0I[1][2] + A1[2][2] * A0I[2][2];
}

inline void ShellElement::getStrainTensor(const double A[3][3],
                                          double E[6]) const
{
    E[0] = 0.5 * (A[0][0] * A[0][0] + A[1][0] * A[1][0] + A[2][0] * A[2][0]) - 0.5;     //E[0][0]
    E[1] = 0.5 * (A[0][1] * A[0][1] + A[1][1] * A[1][1] + A[2][1] * A[2][1]) - 0.5;     //E[1][1]
    E[2] = 0.5 * (A[0][2] * A[0][2] + A[1][2] * A[1][2] + A[2][2] * A[2][2]) - 0.5;     //E[2][2]
    E[3] = 0.5 * (A[0][0] * A[0][1] + A[1][0] * A[1][1] + A[2][0] * A[2][1]);           //E[0][1]
    E[4] = 0.5 * (A[0][0] * A[0][2] + A[1][0] * A[1][2] + A[2][0] * A[2][2]);           //E[0][2]
    E[5] = 0.5 * (A[0][1] * A[0][2] + A[1][1] * A[1][2] + A[2][1] * A[2][2]);           //E[1][2]
}

inline void ShellElement::getFirstDerivatives(double* phi,
                                              double** dphi_dxsi,
                                              const double hm,
                                              const double h0,
                                              const double t,
                                              const double g[3],
                                              const double dt_dxsi[2],
                                              const double A0I[3][3],
                                              const double A[3][3],
                                              const double dg_dxsi[3][2],
                                              double dA_dy[][3][3],
                                              double dE_dy[][6]) const
{
    const int numberOfNodes = base_->getNumberOfNodes();
    for (int i = 0; i < numberOfNodes; ++i)
    {
        double vaux[3];
        vaux[0] = dphi_dxsi[0][i] * A0I[0][0] + dphi_dxsi[1][i] * A0I[1][0];
        vaux[1] = dphi_dxsi[0][i] * A0I[0][1] + dphi_dxsi[1][i] * A0I[1][1];
        vaux[2] = dphi_dxsi[0][i] * A0I[0][2] + dphi_dxsi[1][i] * A0I[1][2];

        // position x
        int ii = 7 * i;
        double (*dA_dy_1)[3] = dA_dy[ii];
        dA_dy_1[0][0] = vaux[0]; dA_dy_1[0][1] = vaux[1]; dA_dy_1[0][2] = vaux[2];
        dA_dy_1[1][0] = 0.0;     dA_dy_1[1][1] = 0.0;     dA_dy_1[1][2] = 0.0;
        dA_dy_1[2][0] = 0.0;     dA_dy_1[2][1] = 0.0;     dA_dy_1[2][2] = 0.0;

        dE_dy[ii][0] =  dA_dy_1[0][0] * A[0][0];
        dE_dy[ii][1] =  dA_dy_1[0][1] * A[0][1];
        dE_dy[ii][2] =  dA_dy_1[0][2] * A[0][2];
        dE_dy[ii][3] = (dA_dy_1[0][0] * A[0][1] + A[0][0] * dA_dy_1[0][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_1[0][0] * A[0][2] + A[0][0] * dA_dy_1[0][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_1[0][1] * A[0][2] + A[0][1] * dA_dy_1[0][2]) * 0.5;

        //position y
        ii = 7 * i + 1;
        double (*dA_dy_2)[3] = dA_dy[ii];
        dA_dy_2[0][0] = 0.0;     dA_dy_2[0][1] = 0.0;     dA_dy_2[0][2] = 0.0;
        dA_dy_2[1][0] = vaux[0]; dA_dy_2[1][1] = vaux[1]; dA_dy_2[1][2] = vaux[2];
        dA_dy_2[2][0] = 0.0;     dA_dy_2[2][1] = 0.0;     dA_dy_2[2][2] = 0.0;

        dE_dy[ii][0] =  dA_dy_2[1][0] * A[1][0];
        dE_dy[ii][1] =  dA_dy_2[1][1] * A[1][1];
        dE_dy[ii][2] =  dA_dy_2[1][2] * A[1][2];
        dE_dy[ii][3] = (dA_dy_2[1][0] * A[1][1] + A[1][0] * dA_dy_2[1][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_2[1][0] * A[1][2] + A[1][0] * dA_dy_2[1][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_2[1][1] * A[1][2] + A[1][1] * dA_dy_2[1][2]) * 0.5;

        // position z
        ii = 7 * i + 2;
        double (*dA_dy_3)[3] = dA_dy[ii];
        dA_dy_3[0][0] = 0.0;     dA_dy_3[0][1] = 0.0;     dA_dy_3[0][2] = 0.0;
        dA_dy_3[1][0] = 0.0;     dA_dy_3[1][1] = 0.0;     dA_dy_3[1][2] = 0.0;
        dA_dy_3[2][0] = vaux[0]; dA_dy_3[2][1] = vaux[1]; dA_dy_3[2][2] = vaux[2];

        dE_dy[ii][0] =  dA_dy_3[2][0] * A[2][0];
        dE_dy[ii][1] =  dA_dy_3[2][1] * A[2][1];
        dE_dy[ii][2] =  dA_dy_3[2][2] * A[2][2];
        dE_dy[ii][3] = (dA_dy_3[2][0] * A[2][1] + A[2][0] * dA_dy_3[2][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_3[2][0] * A[2][2] + A[2][0] * dA_dy_3[2][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_3[2][1] * A[2][2] + A[2][1] * dA_dy_3[2][2]) * 0.5;

        double vaux2[3];
        vaux2[0] = hm * ((1.0 + t * hm) * dphi_dxsi[0][i] + dt_dxsi[0] * phi[i] * hm);
        vaux2[1] = hm * ((1.0 + t * hm) * dphi_dxsi[1][i] + dt_dxsi[1] * phi[i] * hm);
        vaux2[2] = 0.5 * h0 * (1.0 + 2.0 * t * hm) * phi[i];

        vaux[0] = vaux2[0] * A0I[0][0] + vaux2[1] * A0I[1][0] + vaux2[2] * A0I[2][0];
        vaux[1] = vaux2[0] * A0I[0][1] + vaux2[1] * A0I[1][1] + vaux2[2] * A0I[2][1];
        vaux[2] = vaux2[0] * A0I[0][2] + vaux2[1] * A0I[1][2] + vaux2[2] * A0I[2][2];

        // vector x
        ii = 7 * i + 3;
        double (*dA_dy_4)[3] = dA_dy[ii];
        dA_dy_4[0][0] = vaux[0]; dA_dy_4[0][1] = vaux[1]; dA_dy_4[0][2] = vaux[2];
        dA_dy_4[1][0] = 0.0;     dA_dy_4[1][1] = 0.0;     dA_dy_4[1][2] = 0.0;
        dA_dy_4[2][0] = 0.0;     dA_dy_4[2][1] = 0.0;     dA_dy_4[2][2] = 0.0;

        dE_dy[ii][0] =  dA_dy_4[0][0] * A[0][0];
        dE_dy[ii][1] =  dA_dy_4[0][1] * A[0][1];
        dE_dy[ii][2] =  dA_dy_4[0][2] * A[0][2];
        dE_dy[ii][3] = (dA_dy_4[0][0] * A[0][1] + A[0][0] * dA_dy_4[0][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_4[0][0] * A[0][2] + A[0][0] * dA_dy_4[0][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_4[0][1] * A[0][2] + A[0][1] * dA_dy_4[0][2]) * 0.5;

        // vector y
        ii = 7 * i + 4;
        double (*dA_dy_5)[3] = dA_dy[ii];
        dA_dy_5[0][0] = 0.0;     dA_dy_5[0][1] = 0.0;     dA_dy_5[0][2] = 0.0;
        dA_dy_5[1][0] = vaux[0]; dA_dy_5[1][1] = vaux[1]; dA_dy_5[1][2] = vaux[2];
        dA_dy_5[2][0] = 0.0;     dA_dy_5[2][1] = 0.0;     dA_dy_5[2][2] = 0.0;

        dE_dy[ii][0] =  dA_dy_5[1][0] * A[1][0];
        dE_dy[ii][1] =  dA_dy_5[1][1] * A[1][1];
        dE_dy[ii][2] =  dA_dy_5[1][2] * A[1][2];
        dE_dy[ii][3] = (dA_dy_5[1][0] * A[1][1] + A[1][0] * dA_dy_5[1][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_5[1][0] * A[1][2] + A[1][0] * dA_dy_5[1][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_5[1][1] * A[1][2] + A[1][1] * dA_dy_5[1][2]) * 0.5;

        // vector z
        ii = 7 * i + 5;
        double (*dA_dy_6)[3] = dA_dy[ii];
        dA_dy_6[0][0] = 0.0;     dA_dy_6[0][1] = 0.0;     dA_dy_6[0][2] = 0.0;
        dA_dy_6[1][0] = 0.0;     dA_dy_6[1][1] = 0.0;     dA_dy_6[1][2] = 0.0;
        dA_dy_6[2][0] = vaux[0]; dA_dy_6[2][1] = vaux[1]; dA_dy_6[2][2] = vaux[2];

        dE_dy[ii][0] =  dA_dy_6[2][0] * A[2][0];
        dE_dy[ii][1] =  dA_dy_6[2][1] * A[2][1];
        dE_dy[ii][2] =  dA_dy_6[2][2] * A[2][2];
        dE_dy[ii][3] = (dA_dy_6[2][0] * A[2][1] + A[2][0] * dA_dy_6[2][1]) * 0.5;
        dE_dy[ii][4] = (dA_dy_6[2][0] * A[2][2] + A[2][0] * dA_dy_6[2][2]) * 0.5;
        dE_dy[ii][5] = (dA_dy_6[2][1] * A[2][2] + A[2][1] * dA_dy_6[2][2]) * 0.5;

        double maux[3][3];
        maux[0][0] = hm * hm * (phi[i] * dg_dxsi[0][0] + g[0] * dphi_dxsi[0][i]);
        maux[0][1] = hm * hm * (phi[i] * dg_dxsi[0][1] + g[0] * dphi_dxsi[1][i]);
        maux[0][2] = h0 * phi[i] * g[0] * hm;
        maux[1][0] = hm * hm * (phi[i] * dg_dxsi[1][0] + g[1] * dphi_dxsi[0][i]);
        maux[1][1] = hm * hm * (phi[i] * dg_dxsi[1][1] + g[1] * dphi_dxsi[1][i]);
        maux[1][2] = h0 * phi[i] * g[1] * hm;
        maux[2][0] = hm * hm * (phi[i] * dg_dxsi[2][0] + g[2] * dphi_dxsi[0][i]);
        maux[2][1] = hm * hm * (phi[i] * dg_dxsi[2][1] + g[2] * dphi_dxsi[1][i]);
        maux[2][2] = h0 * phi[i] * g[2] * hm;

        // thickness deformation rate
        ii = 7 * i + 6;
        double (*dA_dy_7)[3] = dA_dy[ii];
        dA_dy_7[0][0] = maux[0][0] * A0I[0][0] + maux[0][1] * A0I[1][0] + maux[0][2] * A0I[2][0];
        dA_dy_7[0][1] = maux[0][0] * A0I[0][1] + maux[0][1] * A0I[1][1] + maux[0][2] * A0I[2][1];
        dA_dy_7[0][2] = maux[0][0] * A0I[0][2] + maux[0][1] * A0I[1][2] + maux[0][2] * A0I[2][2];
        dA_dy_7[1][0] = maux[1][0] * A0I[0][0] + maux[1][1] * A0I[1][0] + maux[1][2] * A0I[2][0];
        dA_dy_7[1][1] = maux[1][0] * A0I[0][1] + maux[1][1] * A0I[1][1] + maux[1][2] * A0I[2][1];
        dA_dy_7[1][2] = maux[1][0] * A0I[0][2] + maux[1][1] * A0I[1][2] + maux[1][2] * A0I[2][2];
        dA_dy_7[2][0] = maux[2][0] * A0I[0][0] + maux[2][1] * A0I[1][0] + maux[2][2] * A0I[2][0];
        dA_dy_7[2][1] = maux[2][0] * A0I[0][1] + maux[2][1] * A0I[1][1] + maux[2][2] * A0I[2][1];
        dA_dy_7[2][2] = maux[2][0] * A0I[0][2] + maux[2][1] * A0I[1][2] + maux[2][2] * A0I[2][2];

        dE_dy[ii][0] =  dA_dy_7[0][0] * A[0][0] + dA_dy_7[1][0] * A[1][0] + dA_dy_7[2][0] * A[2][0];
        dE_dy[ii][1] =  dA_dy_7[0][1] * A[0][1] + dA_dy_7[1][1] * A[1][1] + dA_dy_7[2][1] * A[2][1];
        dE_dy[ii][2] =  dA_dy_7[0][2] * A[0][2] + dA_dy_7[1][2] * A[1][2] + dA_dy_7[2][2] * A[2][2];
        dE_dy[ii][3] = (dA_dy_7[0][0] * A[0][1] + dA_dy_7[1][0] * A[1][1] + dA_dy_7[2][0] * A[2][1] +
                        dA_dy_7[0][1] * A[0][0] + dA_dy_7[1][1] * A[1][0] + dA_dy_7[2][1] * A[2][0]) * 0.5;
        dE_dy[ii][4] = (dA_dy_7[0][0] * A[0][2] + dA_dy_7[1][0] * A[1][2] + dA_dy_7[2][0] * A[2][2] +
                        dA_dy_7[0][2] * A[0][0] + dA_dy_7[1][2] * A[1][0] + dA_dy_7[2][2] * A[2][0]) * 0.5;
        dE_dy[ii][5] = (dA_dy_7[0][1] * A[0][2] + dA_dy_7[1][1] * A[1][2] + dA_dy_7[2][1] * A[2][2] +
                        dA_dy_7[0][2] * A[0][1] + dA_dy_7[1][2] * A[1][1] + dA_dy_7[2][2] * A[2][1]) * 0.5;
    }
}

inline void ShellElement::getStrainTensorSecondDerivative(double* phi,
                                                          double** dphi_dxsi,
                                                          const double h0,
                                                          const double hm,
                                                          const int i,
                                                          const int j,
                                                          const double A0I[3][3],
                                                          const double A[3][3],
                                                          const double dA_dy[][3][3],
                                                          double d2E_dydy[6]) const
{
    const int node_a = i / 7;
    const int dof_i = i % 7;
    const int node_b = i / 7;
    const int dof_j = j % 7;
    int index;

    const double (*dA_dy_i)[3] = dA_dy[i];
    const double (*dA_dy_j)[3] = dA_dy[j];

    if ((dof_i >= 3) && (dof_i < 6) && (dof_j == 6))
        index = dof_i - 3;
    else if ((dof_j >= 3) && (dof_j < 6) && (dof_i == 6))
        index = dof_j - 3;
    else
        goto cont1;

    double vaux[3];
    vaux[0] = hm * hm * (phi[node_a] * dphi_dxsi[0][node_b] + phi[node_b] * dphi_dxsi[0][node_a]);
    vaux[1] = hm * hm * (phi[node_a] * dphi_dxsi[1][node_b] + phi[node_b] * dphi_dxsi[1][node_a]);
    vaux[2] = h0 * phi[node_a] * phi[node_b] * hm;

    double d2A_dydy[3];
    d2A_dydy[0] = vaux[0] * A0I[0][0] + vaux[1] * A0I[1][0] + vaux[2] * A0I[2][0];
    d2A_dydy[1] = vaux[0] * A0I[0][1] + vaux[1] * A0I[1][1] + vaux[2] * A0I[2][1];
    d2A_dydy[2] = vaux[0] * A0I[0][2] + vaux[1] * A0I[1][2] + vaux[2] * A0I[2][2];

    d2E_dydy[0] =  dA_dy_i[0][0] * dA_dy_j[0][0] + dA_dy_i[1][0] * dA_dy_j[1][0] + dA_dy_i[2][0] * dA_dy_j[2][0] +
                   d2A_dydy[0] * A[index][0];
    d2E_dydy[1] =  dA_dy_i[0][1] * dA_dy_j[0][1] + dA_dy_i[1][1] * dA_dy_j[1][1] + dA_dy_i[2][1] * dA_dy_j[2][1] +
                  d2A_dydy[1] * A[index][1];
    d2E_dydy[2] =  dA_dy_i[0][2] * dA_dy_j[0][2] + dA_dy_i[1][2] * dA_dy_j[1][2] + dA_dy_i[2][2] * dA_dy_j[2][2] +
                  d2A_dydy[2] * A[index][2];
    d2E_dydy[3] = (dA_dy_i[0][0] * dA_dy_j[0][1] + dA_dy_i[1][0] * dA_dy_j[1][1] + dA_dy_i[2][0] * dA_dy_j[2][1] +
                  dA_dy_j[0][0] * dA_dy_i[0][1] + dA_dy_j[1][0] * dA_dy_i[1][1] + dA_dy_j[2][0] * dA_dy_i[2][1] +
                  d2A_dydy[0] * A[index][1] + A[index][0] * d2A_dydy[1]) * 0.5;
    d2E_dydy[4] = (dA_dy_i[0][0] * dA_dy_j[0][2] + dA_dy_i[1][0] * dA_dy_j[1][2] + dA_dy_i[2][0] * dA_dy_j[2][2] +
                  dA_dy_j[0][0] * dA_dy_i[0][2] + dA_dy_j[1][0] * dA_dy_i[1][2] + dA_dy_j[2][0] * dA_dy_i[2][2] +
                  d2A_dydy[0] * A[index][2] + A[index][0] * d2A_dydy[2]) * 0.5;
    d2E_dydy[5] = (dA_dy_i[0][1] * dA_dy_j[0][2] + dA_dy_i[1][1] * dA_dy_j[1][2] + dA_dy_i[2][1] * dA_dy_j[2][2] +
                  dA_dy_j[0][1] * dA_dy_i[0][2] + dA_dy_j[1][1] * dA_dy_i[1][2] + dA_dy_j[2][1] * dA_dy_i[2][2] +
                  d2A_dydy[1] * A[index][2] + A[index][1] * d2A_dydy[2]) * 0.5;
    return;

    cont1:
    d2E_dydy[0] =  dA_dy_i[0][0] * dA_dy_j[0][0] + dA_dy_i[1][0] * dA_dy_j[1][0] + dA_dy_i[2][0] * dA_dy_j[2][0];
    d2E_dydy[1] =  dA_dy_i[0][1] * dA_dy_j[0][1] + dA_dy_i[1][1] * dA_dy_j[1][1] + dA_dy_i[2][1] * dA_dy_j[2][1];
    d2E_dydy[2] =  dA_dy_i[0][2] * dA_dy_j[0][2] + dA_dy_i[1][2] * dA_dy_j[1][2] + dA_dy_i[2][2] * dA_dy_j[2][2];
    d2E_dydy[3] = (dA_dy_i[0][0] * dA_dy_j[0][1] + dA_dy_i[1][0] * dA_dy_j[1][1] + dA_dy_i[2][0] * dA_dy_j[2][1] +
                  dA_dy_j[0][0] * dA_dy_i[0][1] + dA_dy_j[1][0] * dA_dy_i[1][1] + dA_dy_j[2][0] * dA_dy_i[2][1]) * 0.5;
    d2E_dydy[4] = (dA_dy_i[0][0] * dA_dy_j[0][2] + dA_dy_i[1][0] * dA_dy_j[1][2] + dA_dy_i[2][0] * dA_dy_j[2][2] +
                  dA_dy_j[0][0] * dA_dy_i[0][2] + dA_dy_j[1][0] * dA_dy_i[1][2] + dA_dy_j[2][0] * dA_dy_i[2][2]) * 0.5;
    d2E_dydy[5] = (dA_dy_i[0][1] * dA_dy_j[0][2] + dA_dy_i[1][1] * dA_dy_j[1][2] + dA_dy_i[2][1] * dA_dy_j[2][2] +
                  dA_dy_j[0][1] * dA_dy_i[0][2] + dA_dy_j[1][1] * dA_dy_i[1][2] + dA_dy_j[2][1] * dA_dy_i[2][2]) * 0.5;
}

inline double ShellElement::doubleContraction(const double v1[6],
                                              const double v2[6]) const
{
    return  v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2] +
           (v1[3] * v2[3] + v1[4] * v2[4] + v1[5] * v2[5]) * 2.0;
}