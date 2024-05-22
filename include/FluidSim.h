#ifndef FLUIDSIM_H_
#define FLUIDSIM_H_
#include <cstdlib>
#include <vector>
#include "Particle.h"
#include "FluidCell.h"
#include "Vec3.h"
//#include <string_view>


class FluidSim
{
public:
    // update with changes, such as grid colour
    FluidSim(float _density, float _width, float _height, float _spacing, float _particleRadius, size_t _maxParticles);
    void integrateParticles(float _dt, float _gravity);
    void pushParticlesApart(size_t _numIterations);
    void handleParticleCollisions(float _obstacleX, float _obstacleY, float _obstacleRadius);
    void updateParticleDensity();
    void transferVelocities(bool _toGrid, float _flipRatio);
    void solveIncompressibility(size_t _numIterations, float _dt, float _overRelaxation);
//    void updateParticleColors();
//    void setSciColor(size_t _cellNr, float _val, float _minVal, float _maxVal);
//    void updateCellColors();
    void simulate(float _dt, float _gravity, float _flipRatio, size_t _numIterations, float _overRelaxation, bool _compensateDrift, bool _seperateParticles, float _obstacleX, float _obstacleY, float _obstacleRadius, int step);
    void writeGeo(std::string_view fileName) const;
    void setParticlePos(size_t _index, float _pos);
    void setNumParticles(size_t _numParticles);
    void setS(size_t _index, float _s);
    void setParticles();
    int getfNumX();
    int getfNumY();

private:
    float m_density;
    int m_fNumX;
    int m_fNumY;
    int m_pNumX;
    int m_pNumY;
    float m_h;
    float m_fInvSpacing;
    float m_pInvSpacing;
    int m_fNumCells;
    int m_pNumCells;
    size_t m_numParticles;
    size_t m_maxParticles;

    int const m_fluidCell = 0;
    int const m_airCell = 1;
    int const m_solidCell = 2;

    std::vector<FluidCell> m_fluidCells;
    std::vector<Particle> m_particles;


    // update with the changes we make to the original code
    std::vector<float> m_u;
    std::vector<float> m_v;
    std::vector<float> m_du;
    std::vector<float> m_dv;
    std::vector<float> m_prevU;
    std::vector<float> m_prevV;
    std::vector<float> m_p;
    std::vector<float> m_s;
    std::vector<int> m_cellType;
    std::vector<float> m_cellColor;

    std::vector<float> m_particlePos;
    std::vector<float> m_particleColor;
    std::vector<float> m_particleVel;
    std::vector<float> m_particleDensity;
    float m_particleRestDensity;
    float m_particleRadius;
    std::vector<int> m_numCellParticles;
    std::vector<int> m_firstCellParticle;
    std::vector<int> m_cellParticleIds;
};



#endif