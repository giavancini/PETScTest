#pragma once
#include <vector>

enum class DOFType
{
    POSITION,
    PRESSURE,
    GENERALIZED_VECTOR,
    THICKNESS_RATE
};

class DegreeOfFreedom
{
    public:
        DegreeOfFreedom(const DOFType& type, const double value);

        int getIndex() const;

        DOFType getType() const;

        double getInitialValue() const;

        double getCurrentValue() const;

        double getPastValue() const;

        double getIntermediateValue() const;

        double getInitialFirstTimeDerivative() const;

        double getCurrentFirstTimeDerivative() const;

        double getPastFirstTimeDerivative() const;

        double getIntermediateFirstTimeDerivative() const;

        double getInitialSecondTimeDerivative() const;

        double getCurrentSecondTimeDerivative() const;

        double getPastSecondTimeDerivative() const;

        double getIntermediateSecondTimeDerivative() const;

        void setIndex(const int& index);

        void setInitialValue(const double& value);

        void setCurrentValue(const double& value);

        void setPastValue(const double& value);

        void setIntermediateValue(const double& value);

        void incrementCurrentValue(const double& value);

        void setInitialFirstTimeDerivative(const double& value);

        void setCurrentFirstTimeDerivative(const double& value);

        void setPastFirstTimeDerivative(const double& value);

        void setIntermediateFirstTimeDerivative(const double& value);

        void setInitialSecondTimeDerivative(const double& value);

        void setCurrentSecondTimeDerivative(const double& value);

        void setPastSecondTimeDerivative(const double& value);

        void setIntermediateSecondTimeDerivative(const double& value);

    private:
        int index_;
        DOFType type_;
        double initialValue_;
        double currentValue_;
        double pastValue_;
        double intermediateValue_;
        double initialFirstTimeDerivative_;
        double currentFirstTimeDerivative_;
        double pastFirstTimeDerivative_;
        double intermediateFirstTimeDerivative_;
        double initialSecondTimeDerivative_;
        double currentSecondTimeDerivative_;
        double pastSecondTimeDerivative_;
        double intermediateSecondTimeDerivative_;
};