#include "MathUtil.h"
#include "LogUtil.h"
#include <iostream>
#include <time.h>
// const enum eRotationOrder gRotationOrder = eRotationOrder::XYZ;
// const tVector gGravity = tVector(0, -9.8, 0, 0);
// const tVector gGravity = tVector(0, 0, 0, 0);
cRand cMathUtil::gRand = cRand();

bool cMathUtil::IsPoint(const tVector &vec)
{
    return std::fabs(vec[3] - 1.0) < 1e-10;
}
tVector cMathUtil::VecToPoint(const tVector &vec)
{
    tVector new_vec = vec;
    new_vec[3] = 1;
    return new_vec;
}
int cMathUtil::Clamp(int val, int min, int max)
{
    return std::max(min, std::min(val, max));
}

void cMathUtil::Clamp(const Eigen::VectorXd &min, const Eigen::VectorXd &max,
                      Eigen::VectorXd &out_vec)
{
    out_vec = out_vec.cwiseMin(max).cwiseMax(min);
}

double cMathUtil::Clamp(double val, double min, double max)
{
    return std::max(min, std::min(val, max));
}

double cMathUtil::Saturate(double val) { return Clamp(val, 0.0, 1.0); }

double cMathUtil::Lerp(double t, double val0, double val1)
{
    return (1 - t) * val0 + t * val1;
}

double cMathUtil::NormalizeAngle(double theta)
{
    // normalizes theta to be between [-pi, pi]
    double norm_theta = fmod(theta, 2 * M_PI);
    if (norm_theta > M_PI)
    {
        norm_theta = -2 * M_PI + norm_theta;
    }
    else if (norm_theta < -M_PI)
    {
        norm_theta = 2 * M_PI + norm_theta;
    }
    return norm_theta;
}

double cMathUtil::RandDouble() { return RandDouble(0, 1); }

double cMathUtil::RandDouble(double min, double max)
{
    return gRand.RandDouble(min, max);
}

double cMathUtil::RandDoubleNorm(double mean, double stdev)
{
    return gRand.RandDoubleNorm(mean, stdev);
}

double cMathUtil::RandDoubleExp(double lambda)
{
    return gRand.RandDoubleExp(lambda);
}

double cMathUtil::RandDoubleSeed(double seed)
{
    unsigned int int_seed = *reinterpret_cast<unsigned int *>(&seed);
    std::default_random_engine rand_gen(int_seed);
    std::uniform_real_distribution<double> dist;
    return dist(rand_gen);
}

int cMathUtil::RandInt() { return gRand.RandInt(); }

int cMathUtil::RandInt(int min, int max) { return gRand.RandInt(min, max); }

int cMathUtil::RandUint() { return gRand.RandUint(); }

int cMathUtil::RandUint(unsigned int min, unsigned int max)
{
    return gRand.RandUint(min, max);
}

int cMathUtil::RandIntExclude(int min, int max, int exc)
{
    return gRand.RandIntExclude(min, max, exc);
}

void cMathUtil::SeedRand(unsigned long int seed)
{
    gRand.Seed(seed);
    srand(gRand.RandInt());
}

int cMathUtil::RandSign() { return gRand.RandSign(); }

double cMathUtil::SmoothStep(double t)
{
    double val = t * t * t * (t * (t * 6 - 15) + 10);
    return val;
}

bool cMathUtil::FlipCoin(double p) { return gRand.FlipCoin(p); }


double cMathUtil::Sign(double val) { return SignAux<double>(val); }

int cMathUtil::Sign(int val) { return SignAux<int>(val); }

double cMathUtil::AddAverage(double avg0, int count0, double avg1, int count1)
{
    double total = count0 + count1;
    return (count0 / total) * avg0 + (count1 / total) * avg1;
}

tVector cMathUtil::AddAverage(const tVector &avg0, int count0,
                              const tVector &avg1, int count1)
{
    double total = count0 + count1;
    return (count0 / total) * avg0 + (count1 / total) * avg1;
}

void cMathUtil::AddAverage(const Eigen::VectorXd &avg0, int count0,
                           const Eigen::VectorXd &avg1, int count1,
                           Eigen::VectorXd &out_result)
{
    double total = count0 + count1;
    out_result = (count0 / total) * avg0 + (count1 / total) * avg1;
}

void cMathUtil::CalcSoftmax(const Eigen::VectorXd &vals, double temp,
                            Eigen::VectorXd &out_prob)
{
    assert(out_prob.size() == vals.size());
    int num_vals = static_cast<int>(vals.size());
    double sum = 0;
    double max_val = vals.maxCoeff();
    for (int i = 0; i < num_vals; ++i)
    {
        double val = vals[i];
        val = std::exp((val - max_val) / temp);
        out_prob[i] = val;
        sum += val;
    }

    out_prob /= sum;
}

double cMathUtil::EvalGaussian(const Eigen::VectorXd &mean,
                               const Eigen::VectorXd &covar,
                               const Eigen::VectorXd &sample)
{
    assert(mean.size() == covar.size());
    assert(sample.size() == covar.size());

    Eigen::VectorXd diff = sample - mean;
    double exp_val = diff.dot(diff.cwiseQuotient(covar));
    double likelihood = std::exp(-0.5 * exp_val);

    double partition = CalcGaussianPartition(covar);
    likelihood /= partition;
    return likelihood;
}

double cMathUtil::EvalGaussian(double mean, double covar, double sample)
{
    double diff = sample - mean;
    double exp_val = diff * diff / covar;
    double norm = 1 / std::sqrt(2 * M_PI * covar);
    double likelihood = norm * std::exp(-0.5 * exp_val);
    return likelihood;
}

double cMathUtil::CalcGaussianPartition(const Eigen::VectorXd &covar)
{
    int data_size = static_cast<int>(covar.size());
    double det = covar.prod();
    double partition = std::sqrt(std::pow(2 * M_PI, data_size) * det);
    return partition;
}

double cMathUtil::EvalGaussianLogp(const Eigen::VectorXd &mean,
                                   const Eigen::VectorXd &covar,
                                   const Eigen::VectorXd &sample)
{
    int data_size = static_cast<int>(covar.size());

    Eigen::VectorXd diff = sample - mean;
    double logp = -0.5 * diff.dot(diff.cwiseQuotient(covar));
    double det = covar.prod();
    logp += -0.5 * (data_size * std::log(2 * M_PI) + std::log(det));

    return logp;
}

double cMathUtil::EvalGaussianLogp(double mean, double covar, double sample)
{
    double diff = sample - mean;
    double logp = -0.5 * diff * diff / covar;
    logp += -0.5 * (std::log(2 * M_PI) + std::log(covar));
    return logp;
}

double cMathUtil::Sigmoid(double x) { return Sigmoid(x, 1, 0); }

double cMathUtil::Sigmoid(double x, double gamma, double bias)
{
    double exp = -gamma * (x + bias);
    double val = 1 / (1 + std::exp(exp));
    return val;
}

int cMathUtil::SampleDiscreteProb(const Eigen::VectorXd &probs)
{
    assert(std::abs(probs.sum() - 1) < 0.00001);
    double rand = RandDouble();

    int rand_idx = gInvalidIdx;
    int num_probs = static_cast<int>(probs.size());
    for (int i = 0; i < num_probs; ++i)
    {
        double curr_prob = probs[i];
        rand -= curr_prob;

        if (rand <= 0)
        {
            rand_idx = i;
            break;
        }
    }
    return rand_idx;
}

/**
 * \briewf          categorical random distribution
 */
int cMathUtil::RandIntCategorical(const std::vector<double> &prop)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::discrete_distribution<> d(prop.begin(), prop.end());
    int num = d(gen);
    return num;
}

tVector cMathUtil::CalcBarycentric(const tVector &p, const tVector &a,
                                   const tVector &b, const tVector &c)
{
    tVector v0 = b - a;
    tVector v1 = c - a;
    tVector v2 = p - a;

    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);
    double denom = d00 * d11 - d01 * d01;
    double v = (d11 * d20 - d01 * d21) / denom;
    double w = (d00 * d21 - d01 * d20) / denom;
    double u = 1.0f - v - w;
    return tVector(u, v, w, 0);
}

bool cMathUtil::ContainsAABB(const tVector &pt, const tVector &aabb_min,
                             const tVector &aabb_max)
{
    bool contains = pt[0] >= aabb_min[0] && pt[1] >= aabb_min[1] &&
                    pt[2] >= aabb_min[2] && pt[0] <= aabb_max[0] &&
                    pt[1] <= aabb_max[1] && pt[2] <= aabb_max[2];
    return contains;
}

bool cMathUtil::ContainsAABB(const tVector &aabb_min0, const tVector &aabb_max0,
                             const tVector &aabb_min1, const tVector &aabb_max1)
{
    return ContainsAABB(aabb_min0, aabb_min1, aabb_max1) &&
           ContainsAABB(aabb_max0, aabb_min1, aabb_max1);
}

bool cMathUtil::ContainsAABBXZ(const tVector &pt, const tVector &aabb_min,
                               const tVector &aabb_max)
{
    bool contains = pt[0] >= aabb_min[0] && pt[2] >= aabb_min[2] &&
                    pt[0] <= aabb_max[0] && pt[2] <= aabb_max[2];
    return contains;
}

bool cMathUtil::ContainsAABBXZ(const tVector &aabb_min0,
                               const tVector &aabb_max0,
                               const tVector &aabb_min1,
                               const tVector &aabb_max1)
{
    return ContainsAABBXZ(aabb_min0, aabb_min1, aabb_max1) &&
           ContainsAABBXZ(aabb_max0, aabb_min1, aabb_max1);
}

void cMathUtil::CalcAABBIntersection(const tVector &aabb_min0,
                                     const tVector &aabb_max0,
                                     const tVector &aabb_min1,
                                     const tVector &aabb_max1, tVector &out_min,
                                     tVector &out_max)
{
    out_min = aabb_min0.cwiseMax(aabb_min1);
    out_max = aabb_max0.cwiseMin(aabb_max1);
    if (out_min[0] > out_max[0])
    {
        out_min[0] = 0;
        out_max[0] = 0;
    }
    if (out_min[1] > out_max[1])
    {
        out_min[1] = 0;
        out_max[1] = 0;
    }
    if (out_min[2] > out_max[2])
    {
        out_min[2] = 0;
        out_max[2] = 0;
    }
}

void cMathUtil::CalcAABBUnion(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1, tVector &out_min,
                              tVector &out_max)
{
    out_min = aabb_min0.cwiseMin(aabb_min1);
    out_max = aabb_max0.cwiseMax(aabb_max1);
}

bool cMathUtil::IntersectAABB(const tVector &aabb_min0,
                              const tVector &aabb_max0,
                              const tVector &aabb_min1,
                              const tVector &aabb_max1)
{
    tVector center0 = 0.5 * (aabb_max0 + aabb_min0);
    tVector center1 = 0.5 * (aabb_max1 + aabb_min1);
    tVector size0 = aabb_max0 - aabb_min0;
    tVector size1 = aabb_max1 - aabb_min1;
    tVector test_len = 0.5 * (size0 + size1);
    tVector delta = center1 - center0;
    bool overlap = (std::abs(delta[0]) <= test_len[0]) &&
                   (std::abs(delta[1]) <= test_len[1]) &&
                   (std::abs(delta[2]) <= test_len[2]);
    return overlap;
}

bool cMathUtil::IntersectAABBXZ(const tVector &aabb_min0,
                                const tVector &aabb_max0,
                                const tVector &aabb_min1,
                                const tVector &aabb_max1)
{
    tVector center0 = 0.5 * (aabb_max0 + aabb_min0);
    tVector center1 = 0.5 * (aabb_max1 + aabb_min1);
    tVector size0 = aabb_max0 - aabb_min0;
    tVector size1 = aabb_max1 - aabb_min1;
    tVector test_len = 0.5 * (size0 + size1);
    tVector delta = center1 - center0;
    bool overlap = (std::abs(delta[0]) <= test_len[0]) &&
                   (std::abs(delta[2]) <= test_len[2]);
    return overlap;
}

bool cMathUtil::CheckNextInterval(double delta, double curr_val,
                                  double int_size)
{
    double pad = 0.001 * delta;
    int curr_count = static_cast<int>(std::floor((curr_val + pad) / int_size));
    int prev_count =
        static_cast<int>(std::floor((curr_val + pad - delta) / int_size));
    bool new_action = (curr_count != prev_count);
    return new_action;
}

tVector cMathUtil::SampleRandPt(const tVector &bound_min,
                                const tVector &bound_max)
{
    tVector pt = tVector(RandDouble(bound_min[0], bound_max[0]),
                         RandDouble(bound_min[1], bound_max[1]),
                         RandDouble(bound_min[2], bound_max[2]), 0);
    return pt;
}

tVector cMathUtil::SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max)
{
    return SampleRandPtBias(bound_min, bound_max,
                            0.5 * (bound_max + bound_min));
}

tVector cMathUtil::SampleRandPtBias(const tVector &bound_min,
                                    const tVector &bound_max,
                                    const tVector &focus)
{
    double t = RandDouble(0, 1);
    tVector size = bound_max - bound_min;
    tVector new_min = focus + (t * 0.5) * size;
    tVector new_max = focus - (t * 0.5) * size;
    tVector offset = (bound_min - new_min).cwiseMax(0);
    offset += (bound_max - new_max).cwiseMin(0);
    new_min += offset;
    new_max += offset;

    return SampleRandPt(new_min, new_max);
}

// tQuaternion cMathUtil::RotMatToQuaternion(const tMatrix &mat)
//{
//	//
// http://www.iri.upc.edu/files/scidoc/2068-Accurate-Computation-of-Quaternions-from-Rotation-Matrices.pdf
//	double eta = 0;
//	double q1, q2, q3, q4;	// = [w, x, y, z]
//
//	// determine q1
//	{
//		double detect_value = mat(0, 0) + mat(1, 1) + mat(2, 2);
//		if (detect_value > eta)
//		{
//			q1 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(2, 1) - mat(1, 2), 2);
//			numerator += std::pow(mat(0, 2) - mat(2, 0), 2);
//			numerator += std::pow(mat(1, 0) - mat(0, 1), 2);
//			q1 = 0.5 *  std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q2
//	{
//		double detect_value = mat(0, 0) - mat(1, 1) - mat(2, 2);
//		if (detect_value > eta)
//		{
//			q2 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(2, 1) - mat(1, 2), 2);
//			numerator += std::pow(mat(0, 1) + mat(1, 0), 2);
//			numerator += std::pow(mat(2, 0) + mat(0, 2), 2);
//			q2 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q3
//	{
//		double detect_value = -mat(0, 0) + mat(1, 1) - mat(2, 2);
//		if (detect_value > eta)
//		{
//			q3 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(0, 2) - mat(2, 0), 2);
//			numerator += std::pow(mat(0, 1) + mat(1, 0), 2);
//			numerator += std::pow(mat(1, 2) + mat(2, 1), 2);
//			q3 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	// determine q4
//	{
//		double detect_value = -mat(0, 0) - mat(1, 1) + mat(2, 2);
//		if (detect_value > eta)
//		{
//			q4 = 0.5 * std::sqrt(1 + detect_value);
//		}
//		else
//		{
//			double numerator = 0;
//			numerator += std::pow(mat(1, 0) - mat(0, 1), 2);
//			numerator += std::pow(mat(2, 0) + mat(0, 2), 2);
//			numerator += std::pow(mat(2, 1) + mat(1, 2), 2);
//			q4 = 0.5 * std::sqrt(numerator / (3 - detect_value));
//		}
//	}
//
//	return tQuaternion(q1, q2, q3, q4);
//}
