#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <queue>
#include <chrono>
#include <cstring>
#include <random>

using namespace std;

ifstream fin("input.in");

/// Tim Sort

const int run = 32;

void insertionSort(vector <int> &v, int left, int right)
{
    int i, j;

    if (right >= v.size())
        right = v.size() - 1;

    for (i = left + 1; i <= right; ++i)
    {
        int pos = i;

        for (j = i - 1; j >= left; --j)
            if (v[pos] < v[j])
                swap(v[pos--], v[j]);
    }
}

void merge(vector <int> &v, int left, int mid, int right)
{
    int i = left, j = mid + 1;
    vector <int> aux;

    while (i <= mid && j <= right)
        if (v[i] < v[j]) aux.push_back(v[i++]);
        else aux.push_back(v[j++]);

    while (i <= mid)
        aux.push_back(v[i++]);

    while (j <= right)
        aux.push_back(v[j++]);

    for (i = left, j = 0; i <= right; ++i, ++j)
        v[i] = aux[j];
}

void timSort(vector <int> &v, int n)
{
    int i;

    for (i = 0; i < n; i += run)
        insertionSort(v, i, i + run - 1);

    for (int size = run; size < n; size *= 2)
        for (int left = 0; left + size - 1 < n; left += 2 * size)
        {
            int mid = left + size - 1;
            int right = min(left + 2 * size - 1, n - 1);

            merge(v, left, mid, right);
        }
}

/// Heap Sort

void heapSort(vector <int> &v, int left, int right)
{
    priority_queue <int> heap;

    for (int i = left; i <= right; ++i)
        heap.push(v[i]);

    int i = right;

    while (!heap.empty())
    {
        v[i--] = heap.top();
        heap.pop();
    }
}

/// Intro Sort - Median of Three for Quick Sort

const int RUN = 16;

inline int medianOfThree(vector <int> &v, int left, int mid, int right)
{
    // am folosit xor pentru a determina mediana pentru ca este mai rapid

    if ((v[left] > v[mid]) ^ (v[left] > v[right]))
        return left;
    else
        if ((v[mid] > v[left]) ^ (v[mid] > v[right]))
            return mid;
        else
            return right;
}

int partitionThree(vector <int> &v, int left, int right)
{
    int mid = (left + right) / 2;
    int median = medianOfThree(v, left, mid, right);
    int pivot = v[median];

    if (median != right)
        swap(v[median], v[right]);

    int newPivot = left - 1;

    for (int i = left; i < right; ++i)
        if (v[i] <= pivot)
            swap(v[i], v[++newPivot]);

    swap(v[newPivot+1], v[right]);

    return newPivot + 1;
}

void introThree(vector <int> &v, int left, int right, int depthLimit)
{
    if (right - left + 1 <= RUN)
        insertionSort(v, left, right);
    else
        if (depthLimit == 0)
           heapSort(v, left, right);
        else
        {
            int pivot = partitionThree(v, left, right);

            introThree(v, left, pivot - 1, depthLimit - 1);
            introThree(v, pivot + 1, right, depthLimit - 1);
        }
}

void introSortThree(vector <int> &v, int n)
{
    int depthLimit = 2 * floor(log(n));

    introThree(v, 0, n - 1, depthLimit);
}

/// Intro Sort - Median of Five for Quick Sort

int medianOfFive(int *a, int *b, int *c, int *d, int *e)
{
    int *aux;

    if (*b < *a) { aux = a; a = b; b = aux; }
    if (*d < *c) { aux = c; c = d; d = aux; }
    if (*c < *a) { aux = b; b = d; d = aux; c = a; }

    a = e;

    if (*b < *a) { aux = a; a = b; b = aux; }
    if (*a < *c) { aux = b; b = d; d = aux; a = c; }

    if (*d < *a)
        return *d;
    else
        return *a;
}

int partitionFive(vector <int> &v, int left, int right)
{
    int mid = (left + right) / 2;
    int q1 = (left + mid) / 2, q2 = (mid + right) / 2;
    int pivot = medianOfFive(&v[left], &v[q1], &v[mid], &v[q2], &v[right]);
    int median;

    if (v[left] == pivot) median = left;
    if (v[q1] == pivot) median = q1;
    if (v[mid] == pivot) median = mid;
    if (v[q2] == pivot) median = q2;
    if (v[right] == pivot) median = right;

    if (median != right)
        swap(v[median], v[right]);

    int newPivot = left - 1;

    for (int i = left; i < right; ++i)
        if (v[i] <= pivot)
            swap(v[i], v[++newPivot]);

    swap(v[newPivot+1], v[right]);

    return newPivot + 1;
}

void introFive(vector <int> &v, int left, int right, int depthLimit)
{
    if (right - left + 1 <= run)
        insertionSort(v, left, right);
    else
        if (depthLimit == 0)
            heapSort(v, left, right);
        else
        {
            int pivot = partitionFive(v, left, right);

            introFive(v, left, pivot - 1, depthLimit - 1);
            introFive(v, pivot + 1, right, depthLimit - 1);
        }
}

void introSortFive(vector <int> &v, int n)
{
    int depthLimit = 2 * floor(log(n));

    introFive(v, 0, n - 1, depthLimit);
}

/// Radix Sort 2^4

void radixSort4(vector <int> &v, int n)
{
    // de i ori shiftam cu 4 biti la dreapta
    for (int i = 0; i < 8; ++i)
    {
        vector <int> bucket[16];

        for (int j = 0; j < n; ++j)
        {
            int digit = v[j] >> (4 * i) & 15;
            bucket[digit].push_back(v[j]);
        }

        int index = 0;

        for (int j = 0; j < 16; ++j)
            for (int k = 0; k < bucket[j].size(); ++k)
                v[index++] = bucket[j][k];
    }
}

/// Radix Sort 2^8

void radixSort8(vector <int> &v, int n)
{
    // de i ori shiftam cu 8 biti la dreapta
    for (int i = 0; i < 4; ++i)
    {
        vector <int> bucket[256];

        for (int j = 0; j < n; ++j)
        {
            int digit = v[j] >> (8 * i) & 255;
            bucket[digit].push_back(v[j]);
        }

        int index = 0;

        for (int j = 0; j < 256; ++j)
            for (int k = 0; k < bucket[j].size(); ++k)
                v[index++] = bucket[j][k];
    }
}

/// other functions and main

int strToNum(char *p)
{
    int num = 0;

    for (int i = 0; p[i]; ++i)
        num = num * 10 + p[i] - '0';

    return num;
}

void generateVector(vector <int> &v, int n, int max)
{
    random_device device;
    mt19937 generator(device());
    uniform_int_distribution <int> distribution(0,max);

    for (int i = 0; i < n; ++i)
        v[i] = distribution(generator);
}

void copyVector(vector <int> &v, vector <int> &aux, int n)
{
    for (int i = 0; i < n; ++i)
        aux[i] = v[i];
}

bool isSorted(vector <int> &v)
{
    for (int i = 1; i < v.size(); ++i)
        if (v[i] < v[i-1])
            return false;

    return true;
}

void callSort(vector <int> &v, int n, char *s)
{
    if (!strcmp(s, "insertion"))
    {
        if (n >= 100000)
            cout << "Can't perform Insertion Sort! Data set is too large.\n";
        else
        {
            auto startTime = std::chrono::high_resolution_clock::now();
            insertionSort(v, 0, n);
            auto endTime = std::chrono::high_resolution_clock::now();
            auto time = endTime - startTime;

            if (isSorted(v))
                cout << "SUCCESS! Insertion Sort took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
            else
                cout << "OH NO! Insertion Sort failed :(\n";
        }
    }

    if (!strcmp(s, "heap"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        heapSort(v, 0, n - 1);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Heap Sort took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Heap Sort failed :(\n";
    }

    if (!strcmp(s, "radix4"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        radixSort4(v, n);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Radix Sort (base 2^4) took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Radix Sort (base 2^4) failed :(\n";
    }

    if (!strcmp(s, "radix8"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        radixSort8(v, n);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Radix Sort (base 2^8) took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Radix Sort (base 2^8) failed :(\n";
    }

    if (!strcmp(s, "intro3"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        introSortThree(v, n);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Intro Sort (with the median of three as pivot for Quick Sort) took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Intro Sort (with the median of three as pivot for Quick Sort) failed :(\n";
    }

    if (!strcmp(s, "intro5"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        introSortFive(v, n);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Intro Sort (with the median of five as pivot for Quick Sort) took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Intro Sort (with the median of five as pivot for Quick Sort) failed :(\n";
    }

    if (!strcmp(s, "tim"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        timSort(v, n);
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! Tim Sort took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            cout << "OH NO! Tim Sort failed :(\n";
    }

    if (!strcmp(s, "STL"))
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        sort(v.begin(), v.end());
        auto endTime = std::chrono::high_resolution_clock::now();
        auto time = endTime - startTime;

        if (isSorted(v))
            cout << "SUCCESS! STL Sort took " << time / std::chrono::milliseconds(1) / 1000.0 << " seconds.\n";
        else
            // unlikely, but
            cout << "OH NO! STL Sort failed :(\n";
    }
}

int main()
{
    char s[100], *p;

    fin.getline(s, 100);

    p = strtok(s, " ");
    p = strtok(NULL, " ");
    p = strtok(NULL, " ");

    int tests = strToNum(p);
    int n, max;
    char sorts[8][20];

    strcpy(sorts[0], "insertion");
    strcpy(sorts[1], "heap");
    strcpy(sorts[2], "radix4");
    strcpy(sorts[3], "radix8");
    strcpy(sorts[4], "intro3");
    strcpy(sorts[5], "intro5");
    strcpy(sorts[6], "tim");
    strcpy(sorts[7], "STL");

    for (int i = 0; i < tests; ++i)
    {
        fin.getline(s, 100);

        p = strtok(s, " ");
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");

        n = strToNum(p);

        p = strtok(NULL, " ");
        p = strtok(NULL, " ");
        p = strtok(NULL, " ");

        max = strToNum(p);

        cout << "Test " << i + 1 << '\n' << "N = " << n << " Max = " << max << '\n';

        vector <int> v(n);
        vector <int> aux(n);

        generateVector(v, n, max);
        copyVector(v, aux, n);

        for (int j = 0; j < 8; ++j)
        {
            copyVector(aux, v, n);
            callSort(v, n, sorts[j]);
        }

        cout << '\n';
    }

    return 0;
}