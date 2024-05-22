#include <iostream>
#include <cmath>
#include <algorithm>
#include <fmt/format.h>
#include <fstream>
#include "FluidSim.h"

FluidSim::FluidSim(float _density, float _width, float _height, float _spacing, float _particleRadius, size_t _maxParticles)
{
    m_density = _density;
    m_fNumX = std::floor(_width / _spacing) + 1;
    m_fNumY = std::floor(_height / _spacing) + 1;
    m_h = std::max(_width / m_fNumX, _height / m_fNumY);
    m_fInvSpacing = 1.0f / m_h;
    m_fNumCells = m_fNumX * m_fNumY;

    m_u.resize(m_fNumCells);
    m_v.resize(m_fNumCells);
    m_du.resize(m_fNumCells);
    m_dv.resize(m_fNumCells);
    m_prevU.resize(m_fNumCells);
    m_prevV.resize(m_fNumCells);
    m_p.resize(m_fNumCells);
    m_s.resize(m_fNumCells);
    m_cellType.resize(m_fNumCells);
    m_cellColor.resize(3 * m_fNumCells);

    m_maxParticles = _maxParticles;

    m_particlePos.resize(2 * m_maxParticles);
    m_particleColor.resize(3 * m_maxParticles);
    for (size_t i=0; i<m_maxParticles; i++)
    {
        m_particleColor[3 * i + 2] = 1.0;
    }
    m_particleVel.resize(2*m_maxParticles);
    m_particleDensity.resize(m_fNumCells);
    m_particleRestDensity = 0.0;

    m_particleRadius = _particleRadius;
    m_pInvSpacing = 1.0f/ (2.2f * _particleRadius);
    m_pNumX = std::floor(_width * m_pInvSpacing) + 1;
    m_pNumY = std::floor(_height / m_pInvSpacing) + 1;
    m_pNumCells = m_pNumX * m_pNumY;

    m_numCellParticles.resize(m_pNumCells);
    m_firstCellParticle.resize(m_pNumCells + 1);
    m_cellParticleIds.resize(m_maxParticles);

    m_numParticles = 0;
}

void FluidSim::integrateParticles(float _dt, float _gravity)
{
    for (size_t i=0; i<m_numParticles; i++)
    {
        m_particleVel[2 * i + 1] += _dt * _gravity;
        m_particlePos[2 * i] += m_particleVel[2 * i] * _dt;
        m_particlePos[2 * i + 1] += m_particleVel[2 * i + 1] * _dt;

//        m_particles[i].vel.y += _dt * _gravity;
//        m_particles[i].pos.x =
    }
}

void FluidSim::pushParticlesApart(size_t _numIterations)
{
    float colorDiffusionCoeff = 0.001;

    // count particles in each cell

    m_numCellParticles.assign(m_pNumCells, 0);
    for (size_t i=0; i<m_numParticles; i++)
    {
        float x = m_particlePos[2 * i];
        float y = m_particlePos[2 * i + 1];

        int xi = std::clamp(static_cast<int>(std::floor(x * m_pInvSpacing)), 0, static_cast<int>(m_pNumX - 1));
        int yi = std::clamp(static_cast<int>(std::floor(y * m_pInvSpacing)), 0, static_cast<int>(m_pNumY - 1));

        int cellNr = xi * static_cast<int>(m_pNumY) + yi;
        m_firstCellParticle[cellNr]--;
        m_cellParticleIds[m_firstCellParticle[cellNr]] = i;
    }

    // push particles apart

    float minDist = 2.0 * m_particleRadius;
    float minDist2 = minDist * minDist;

    for (size_t iter = 0; iter < _numIterations; iter++)
    {
        for (size_t i=0; i<m_numParticles; i++)
        {
            float px = m_particlePos[2 * i];
            float py = m_particlePos[2 * i + 1];

            int pxi = static_cast<int>(std::floor(px * m_pInvSpacing));
            int pyi = static_cast<int>(std::floor(py * m_pInvSpacing));
            int x0 = std::max(static_cast<int>(pxi - 1), 0);
            int y0 = std::max(static_cast<int>(pyi - 1), 0);
            int x1 = std::min(static_cast<int>(pxi + 1), static_cast<int>(m_pNumX - 1));
            int y1 = std::min(static_cast<int>(pyi + 1), static_cast<int>(m_pNumY - 1));

            for (int xi = x0; xi <= x1; xi++) {
                for (int yi = y0; yi <= y1; yi++) {
                    int cellNr = xi * m_pNumY + yi;
                    int first = m_firstCellParticle[cellNr];
                    int last = m_firstCellParticle[cellNr + 1];
                    for (int j = first; j < last; j++) {
                        int id = m_cellParticleIds[j];
                        if (id == i)
                            continue;
                        float qx = m_particlePos[2 * id];
                        float qy = m_particlePos[2 * id + 1];

                        float dx = qx - px;
                        float dy = qy - py;
                        float d2 = dx * dx + dy * dy;
                        if (d2 > minDist2 || d2 == 0.0f)
                            continue;
                        float d = std::sqrt(d2);
                        float s = 0.5f * (minDist - d) / d;
                        dx *= s;
                        dy *= s;
                        m_particlePos[2 * i] -= dx;
                        m_particlePos[2 * i + 1] -= dy;
                        m_particlePos[2 * id] += dx;
                        m_particlePos[2 * id + 1] += dy;

                        // Diffuse colors

                        for (int k = 0; k < 3; k++) {
                            float color0 = m_particleColor[3 * i + k];
                            float color1 = m_particleColor[3 * id + k];
                            float color = (color0 + color1) * 0.5f;
                            m_particleColor[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
                            m_particleColor[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
                        }
                    }
                }
            }

        }
    }
}

void FluidSim::handleParticleCollisions(float _obstacleX, float _obstacleY, float _obstacleRadius)
{
    float h = 1.0f / m_fInvSpacing;
    float r = m_particleRadius;
    float or_ = _obstacleRadius;
    float or2 = or_ * or_;
    float minDist = _obstacleRadius + r;
    float minDist2 = minDist * minDist;

    float minX = h + r;
    float maxX = (m_fNumX - 1) * h - r;
    float minY = h + r;
    float maxY = (m_fNumY - 1) * h - r;

    for (int i = 0; i < m_numParticles; i++) {
        float x = m_particlePos[2 * i];
        float y = m_particlePos[2 * i + 1];

        float dx = x - _obstacleX;
        float dy = y - _obstacleY;
        float d2 = dx * dx + dy * dy;

        // Obstacle collision

//        if (d2 < minDist2) {
//            // Handle obstacle collision
//            m_particleVel[2 * i] = scene.obstacleVelX;
//            m_particleVel[2 * i + 1] = scene.obstacleVelY;
//        }

        // Wall collisions

        if (x < minX) {
            x = minX;
            m_particleVel[2 * i] = 0.0f;
        }
        if (x > maxX) {
            x = maxX;
            m_particleVel[2 * i] = 0.0f;
        }
        if (y < minY) {
            y = minY;
            m_particleVel[2 * i + 1] = 0.0f;
        }
        if (y > maxY) {
            y = maxY;
            m_particleVel[2 * i + 1] = 0.0f;
        }

        // Update particle positions
        m_particlePos[2 * i] = x;
        m_particlePos[2 * i + 1] = y;
    }

}

void FluidSim::updateParticleDensity()
{
    int n = m_fNumY;
    float h = m_h;
    float h1 = m_fInvSpacing;
    float h2 = 0.5f * h;

    auto& d = m_particleDensity;

    std::fill(d.begin(), d.end(), 0.0f);

    for (int i = 0; i < m_numParticles; i++) {
        float x = m_particlePos[2 * i];
        float y = m_particlePos[2 * i + 1];

        x = std::clamp(x, h, (m_fNumX - 1) * h);
        y = std::clamp(y, h, (m_fNumY - 1) * h);

        int x0 = static_cast<int>(std::floor((x - h2) * h1));
        float tx = ((x - h2) - x0 * h) * h1;
        int x1 = std::min(x0 + 1, static_cast<int>(m_fNumX - 2));

        int y0 = static_cast<int>(std::floor((y - h2) * h1));
        float ty = ((y - h2) - y0 * h) * h1;
        int y1 = std::min(y0 + 1, static_cast<int>(m_fNumY - 2));

        float sx = 1.0f - tx;
        float sy = 1.0f - ty;

        if (x0 < m_fNumX && y0 < m_fNumY) d[x0 * n + y0] += sx * sy;
        if (x1 < m_fNumX && y0 < m_fNumY) d[x1 * n + y0] += tx * sy;
        if (x1 < m_fNumX && y1 < m_fNumY) d[x1 * n + y1] += tx * ty;
        if (x0 < m_fNumX && y1 < m_fNumY) d[x0 * n + y1] += sx * ty;
    }

    if (m_particleRestDensity == 0.0f) {
        float sum = 0.0f;
        int numFluidCells = 0;

        for (int i = 0; i < m_fNumCells; i++) {
            if (m_cellType[i] == m_fluidCell) {
                sum += d[i];
                numFluidCells++;
            }
        }

        if (numFluidCells > 0)
            m_particleRestDensity = sum / numFluidCells;
    }
}

void FluidSim::transferVelocities(bool _toGrid, float _flipRatio)
{
    int n = m_fNumY;
    float h = m_h;
    float h1 = m_fInvSpacing;
    float h2 = 0.5f * h;

    if (_toGrid)
    {
        m_prevU = m_u;
        m_prevV = m_v;

        std::fill(m_du.begin(), m_du.end(), 0.0f);
        std::fill(m_dv.begin(), m_dv.end(), 0.0f);
        std::fill(m_u.begin(), m_u.end(), 0.0f);
        std::fill(m_v.begin(), m_v.end(), 0.0f);

        for (int i = 0; i < m_fNumCells; i++)
            m_cellType[i] = m_s[i] == 0.0f ? m_solidCell : m_airCell;

        for (int i = 0; i < m_numParticles; i++)
        {
            float x = m_particlePos[2 * i];
            float y = m_particlePos[2 * i + 1];
            int xi = std::clamp(static_cast<int>(std::floor(x * h1)), 0, static_cast<int>(m_fNumX) - 1);
            int yi = std::clamp(static_cast<int>(std::floor(y * h1)), 0, static_cast<int>(m_fNumY) - 1);
            int cellNr = xi * n + yi;
            if (m_cellType[cellNr] == m_airCell)
                m_cellType[cellNr] = m_fluidCell;
        }
    }

    for (int component = 0; component < 2; component++)
    {
        float dx = component == 0 ? 0.0f : h2;
        float dy = component == 0 ? h2 : 0.0f;

        std::vector<float> &f = component == 0 ? m_u : m_v;
        std::vector<float> &prevF = component == 0 ? m_prevU : m_prevV;
        std::vector<float> &d = component == 0 ? m_du : m_dv;

        for (int i = 0; i < m_numParticles; i++)
        {
            float x = m_particlePos[2 * i];
            float y = m_particlePos[2 * i + 1];

            x = std::clamp(x, h, (m_fNumX - 1) * h);
            y = std::clamp(y, h, (m_fNumY - 1) * h);

            int x0 = std::min(static_cast<int>(std::floor((x - dx) * h1)), static_cast<int>(m_fNumX - 2));
            float tx = ((x - dx) - x0 * h) * h1;
            int x1 = std::min(x0 + 1, static_cast<int>(m_fNumX - 2));

            int y0 = std::min(static_cast<int>(std::floor((y - dy) * h1)), static_cast<int>(m_fNumY - 2));
            float ty = ((y - dy) - y0 * h) * h1;
            int y1 = std::min(y0 + 1, static_cast<int>(m_fNumY - 2));

            float sx = 1.0f - tx;
            float sy = 1.0f - ty;

            float d0 = sx * sy;
            float d1 = tx * sy;
            float d2 = tx * ty;
            float d3 = sx * ty;

            int nr0 = x0 * n + y0;
            int nr1 = x1 * n + y0;
            int nr2 = x1 * n + y1;
            int nr3 = x0 * n + y1;

            if (_toGrid)
            {
                float pv = m_particleVel[2 * i + component];
                f[nr0] += pv * d0;
                d[nr0] += d0;
                f[nr1] += pv * d1;
                d[nr1] += d1;
                f[nr2] += pv * d2;
                d[nr2] += d2;
                f[nr3] += pv * d3;
                d[nr3] += d3;
            } else
            {
                int offset = component == 0 ? n : 1;
                float valid0 = (m_cellType[nr0] != m_airCell || m_cellType[nr0 - offset] != m_airCell) ? 1.0f : 0.0f;
                float valid1 = (m_cellType[nr1] != m_airCell || m_cellType[nr1 - offset] != m_airCell) ? 1.0f : 0.0f;
                float valid2 = (m_cellType[nr2] != m_airCell || m_cellType[nr2 - offset] != m_airCell) ? 1.0f : 0.0f;
                float valid3 = (m_cellType[nr3] != m_airCell || m_cellType[nr3 - offset] != m_airCell) ? 1.0f : 0.0f;

                float v = m_particleVel[2 * i + component];
                float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

                if (d > 0.0f)
                {
                    float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] +
                                  valid3 * d3 * f[nr3]) / d;
                    float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) +
                                  valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
                    float flipV = v + corr;

                    m_particleVel[2 * i + component] = (1.0f - _flipRatio) * picV + _flipRatio * flipV;
                }
            }
        }

        if (_toGrid)
        {
            for (size_t i = 0; i < f.size(); i++)
            {
                if (d[i] > 0.0f)
                    f[i] /= d[i];
            }

            for (int i = 0; i < m_fNumX; i++)
            {
                for (int j = 0; j < m_fNumY; j++)
                {
                    bool solid = m_cellType[i * n + j] == m_solidCell;
                    if (solid || (i > 0 && m_cellType[(i - 1) * n + j] == m_solidCell))
                        m_u[i * n + j] = m_prevU[i * n + j];
                    if (solid || (j > 0 && m_cellType[i * n + j - 1] == m_solidCell))
                        m_v[i * n + j] = m_prevV[i * n + j];
                }
            }
        }
    }
}

void FluidSim::solveIncompressibility(size_t _numIterations, float _dt, float _overRelaxation)
{
    m_p.assign(m_fNumCells, 0.0f);
    m_prevU = m_u;
    m_prevV = m_v;

    int n = m_fNumY;
    float cp = m_density * m_h / _dt;

    for (int i = 0; i < m_fNumCells; i++) {
        float u = m_u[i];
        float v = m_v[i];
    }

    for (int iter = 0; iter < _numIterations; iter++) {
        for (int i = 1; i < m_fNumX - 1; i++) {
            for (int j = 1; j < m_fNumY - 1; j++) {
                if (m_cellType[i * n + j] != m_fluidCell)
                    continue;

                int center = i * n + j;
                int left = (i - 1) * n + j;
                int right = (i + 1) * n + j;
                int bottom = i * n + j - 1;
                int top = i * n + j + 1;

                float s = m_s[center];
                float sx0 = m_s[left];
                float sx1 = m_s[right];
                float sy0 = m_s[bottom];
                float sy1 = m_s[top];
                s = sx0 + sx1 + sy0 + sy1;
                if (s == 0.0f)
                    continue;

                float div = m_u[right] - m_u[center] + m_v[top] - m_v[center];

                if (m_particleRestDensity > 0.0f) {
                    float k = 1.0f;
                    float compression = m_particleDensity[i * n + j] - m_particleRestDensity;
                    if (compression > 0.0f)
                        div = div - k * compression;
                }

                float p = -div / s;
                p *= _overRelaxation;
                m_p[center] += cp * p;

                m_u[center] -= sx0 * p;
                m_u[right] += sx1 * p;
                m_v[center] -= sy0 * p;
                m_v[top] += sy1 * p;
            }
        }
    }
}

void FluidSim::simulate(float _dt, float _gravity, float _flipRatio, size_t _numIterations, float _overRelaxation,
                        bool _compensateDrift, bool _seperateParticles, float _obstacleX, float _obstacleY,
                        float _obstacleRadius, int step)
{
    int subSteps = 1.0;
    float sdt = _dt / static_cast<float>(subSteps);

//    for (int step=0; step<subSteps; step++)
//    {
        integrateParticles(sdt, _gravity);
        pushParticlesApart(_numIterations);
        handleParticleCollisions(_obstacleX, _obstacleY, _obstacleRadius);

        // in original only had the one argument
        transferVelocities(true, 1);

        updateParticleDensity();
        solveIncompressibility(_numIterations, sdt, _overRelaxation);
        transferVelocities(false, _flipRatio);

//        writeGeo(fmt::format("/transfer/fluid/particle.{:04d}.geo", step));
//    }
}

void FluidSim::writeGeo(std::string_view fileName) const
{
    auto file=std::ofstream(fileName.data());
    file<<"PGEOMETRY V5\n";
    file<<"NPoints "<<m_numParticles<<" NPrims 1\n";
    file<<"NPointGroups 0 NPrimGroups 0\n";
    file<<"NPointAttrib 2 NVertexAttrib 0 NPrimAttrib 1 NAttrib 0\n";
    file<<"PointAttrib\n";
    file<<"Cd 3 float 1 1 1\n";
    file<<"pscale 1 float 0.5\n";
    size_t numParts=0;
//    for (auto p:m_particles)
    for (int i=0; i<m_numParticles; i++)
    {
//        if (p.alive == ParticleState::Alive)
//        {
//            file << p.pos.x << " " << p.pos.y << " " << p.pos.z << " 1.0 (";
//            file << p.colour.x << " " << p.colour.y << " " << p.colour.z << " ";
//            file << p.size << ")\n";
//            ++numParts;

//        }
        file << m_particlePos[2 * i] << " " << m_particlePos[2 * i + 1] << " " << 0 << " 1.0 (";
        file << 0 << " " << 0 << " " << 0 << " ";
        file << 1.0 << ")\n";
    }

    file<<"PrimitiveAttrib\n";
    file<<"generator 1 index 1 papi\n";
    file<<"Part "<<m_numParticles<<" ";
    for(size_t i=0; i<m_numParticles; ++i)
    {
        file<<i<<" ";
    }
    file<<"[0]\n";
    file<<"beginExtra\nendExtra\n";

//    std::cout<< fileName.data() << "\n";

    file.close();
}

void FluidSim::setNumParticles(size_t _numParticles)
{
    m_numParticles = _numParticles;
}

void FluidSim::setParticlePos(size_t _index, float _pos)
{
    m_particlePos[_index] = _pos;
}

void FluidSim::setS(size_t _index, float _s)
{
    m_s[_index] = _s;
}

void FluidSim::setParticles()
{
    m_particles.resize(m_numParticles);
}

int FluidSim::getfNumX()
{
    return m_fNumX;
}

int FluidSim::getfNumY()
{
    return m_fNumY;
}