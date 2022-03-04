#pragma once

class AnalysisParameters
{
    public:
    AnalysisParameters();

    ~AnalysisParameters();

    void setDimension(const int dimension);

    void setNumberOfSteps(const int numberOfSteps);

    void setMaxNonlinearIterations(const int maxNonlinearIterations);

    void setNonlinearTolerance(const double nonlinearTolerance);

    void setDeltat(const double deltat);

    void setCurrentTime(const double time);

    void setGravity(const double gravity_x,
                    const double gravity_y,
                    const double gravity_z = 0.0);
    
    void setSpectralRadius(const double rhoInf);

    void setGeneralizedAlphas(const double alphaM, const double alphaF);

    void setAlpha(const double alpha);

    void setMeshLength(const double h);

    void setModelVolume(const double volume);

    void setInitialAccel(const bool initialAccel);

    void setExportFrequency(const int& freq);

    void setStaticAnalysis(const bool& isStaticAnalysis);

    void setLumpedMass(const bool& useLumpedMass);

    int getDimension() const;
    
    int getNumberOfSteps() const;

    int getMaxNonlinearIterations() const;

    double getNonlinearTolerance() const;

    double getDeltat() const;

    double getCurrentTime() const;

    double* getGravity();

    double getSpectralRadius() const;

    double getAlphaM() const;

    double getAlphaF() const;

    double getGamma() const;

    double getBeta() const;

    double getAlpha() const;

    double getMeshLength() const;

    double getModelVolume() const;

    int getExportFrequency() const;

    bool getInitialAccel() const;

    bool isStaticAnalysis() const;

    bool useLumpedMass() const;

    private:
    int dimension_;
    int numberOfSteps_;
	int maxNonlinearIterations_;
	double nonlinearTolerance_;
	double deltat_;
    double currentTime_;
	double gravity_[3];
	double rhoInf_;
    double alphaM_;
    double alphaF_;
    double gamma_;
    double beta_;
	double alpha_;
	double meshLength_;
    double modelVolume_;
    int exportFrequency_;
	bool initialAccel_;
    bool isStaticAnalysis_;
    bool useLumpedMass_;
};