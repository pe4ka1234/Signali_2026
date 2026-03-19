#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <set>
#include <sstream>
#include <iomanip>
#include <numeric>

using Vector = std::vector<int>;
using Matrix = std::vector<Vector>;

struct DualCodeword
{
    Vector coeffs;
    Vector word;
    int w;
};

Matrix buildHammingParityCheckMatrix(int m)
{
    int n = (1 << m) - 1;
    Matrix H(m, Vector(n, 0));

    for (int x = 1; x <= n; ++x)
    {
        for (int i = 0; i < m; ++i)
        {
            H[i][x - 1] = (x >> i) & 1;
        }
    }

    return H;
}

Vector xorVectors(const Vector& a, const Vector& b)
{
    Vector result(a.size(), 0);
    for (size_t i = 0; i < a.size(); ++i)
    {
        result[i] = a[i] ^ b[i];
    }
    return result;
}

std::string vectorToString(const Vector& v)
{
    std::string s;
    s.reserve(v.size());
    for (int bit : v)
    {
        s += static_cast<char>('0' + bit);
    }
    return s;
}

int weight(const Vector& v)
{
    return std::accumulate(v.begin(), v.end(), 0);
}

std::vector<DualCodeword> generateDualCodewords(const Matrix& H)
{
    int m = static_cast<int>(H.size());
    int n = static_cast<int>(H[0].size());
    int total = 1 << m;

    std::vector<DualCodeword> result;
    result.reserve(total);

    for (int mask = 0; mask < total; ++mask)
    {
        Vector coeffs(m, 0);
        Vector word(n, 0);

        for (int i = 0; i < m; ++i)
        {
            coeffs[i] = (mask >> (m - 1 - i)) & 1;
        }

        for (int i = 0; i < m; ++i)
        {
            if (coeffs[i] == 1)
            {
                word = xorVectors(word, H[i]);
            }
        }

        result.push_back({ coeffs, word, weight(word) });
    }

    return result;
}

std::string buildReport(int m)
{
    Matrix H = buildHammingParityCheckMatrix(m);
    std::vector<DualCodeword> codewords = generateDualCodewords(H);

    int n = (1 << m) - 1;
    int expectedWeight = 1 << (m - 1);

    std::ostringstream out;

    out << "CHECK: dual code of Hamming code for m = " << m << "\n";
    out << "Code length n = 2^" << m << " - 1 = " << n << "\n";
    out << "Expected nonzero weight for simplex code: " << expectedWeight << "\n\n";

    out << "Parity-check matrix H:\n";
    for (const auto& row : H)
    {
        out << "  " << vectorToString(row) << "\n";
    }

    out << "\nDual codewords:\n\n";

    std::set<int> nonzeroWeights;

    for (size_t i = 0; i < codewords.size(); ++i)
    {
        const auto& cw = codewords[i];

        if (cw.w != 0)
        {
            nonzeroWeights.insert(cw.w);
        }

        out << std::setw(2) << (i + 1)
            << ") coeffs: " << vectorToString(cw.coeffs)
            << "   word: " << vectorToString(cw.word)
            << "   weight: " << cw.w << "\n";
    }

    out << "\nSummary:\n";
    out << "  Total number of codewords: " << codewords.size() << "\n";
    out << "  Unique nonzero weights: ";

    if (nonzeroWeights.empty())
    {
        out << "none";
    }
    else
    {
        out << "[";
        bool first = true;
        for (int w : nonzeroWeights)
        {
            if (!first) out << ", ";
            out << w;
            first = false;
        }
        out << "]";
    }
    out << "\n";

    if (nonzeroWeights.size() == 1 && *nonzeroWeights.begin() == expectedWeight)
    {
        out << "  Conclusion: all nonzero codewords have weight " << expectedWeight << ".\n";
        out << "  Therefore, the dual code of the Hamming code is indeed a simplex code.\n";
    }
    else
    {
        out << "  Conclusion: the simplex code property is NOT satisfied.\n";
    }

    return out.str();
}

void saveToFile(const std::string& filename, const std::string& text)
{
    std::ofstream file(filename);
    if (!file)
    {
        std::cerr << "Error: could not open file " << filename << " for writing.\n";
        return;
    }

    file << text;
    file.close();
}

int main()
{
    int m;

    std::cout << "Enter m for Hamming code ";
    std::cin >> m;

    if (std::cin.fail() || m < 2 || m > 10)
    {
        std::cout << "Incorrect input. Enter an integer m from 2 to 10.\n";
        return 1;
    }

    std::string report = buildReport(m);

    std::cout << "\n" << report << "\n";

    std::string filename = "dual_hamming_m" + std::to_string(m) + ".txt";
    saveToFile(filename, report);

    std::cout << "\nThe result was also saved to file: " << filename << "\n";

    return 0;
}