#include <iostream>
#include <vector>
#include <chrono>
#include <thread>
#include <cmath>

//COmplex struct - can handle: addition, multiplication, division, inversion,
//absolute value, double conversion, adaptive printing for pure real and imaginary aswell
//NOT ABLE TO HANDLE - euler represenatation - direct argument calculation 
struct complex {double RealPart;double ImaginaryPart;
    complex(double r = 0.0, double i = 0.0) : RealPart(r), ImaginaryPart(i) {}
//complex addition
    complex operator+(const complex& other) const {
        return {RealPart + other.RealPart, ImaginaryPart + other.ImaginaryPart};
    }
    complex operator+(double val) const { return {RealPart + val, ImaginaryPart}; }
    friend complex operator+(double val, const complex& c) { return c + val; }

    // --- Subtraction ---
    complex operator-(const complex& other) const {
        return {RealPart - other.RealPart, ImaginaryPart - other.ImaginaryPart};}
//complex subtraction
    complex operator-(double val) const { return {RealPart - val, ImaginaryPart}; }
    friend complex operator-(double val, const complex& c) { return complex(val) - c; }
//complex multiplication
    complex operator*(const complex& other) const {
        return {RealPart * other.RealPart - ImaginaryPart * other.ImaginaryPart,
                RealPart * other.ImaginaryPart + ImaginaryPart * other.RealPart};
    }
    complex operator*(double val) const { return {RealPart * val, ImaginaryPart * val}; }
    friend complex operator*(double val, const complex& c) { return c * val; }
//complex division and inversion
    complex Invert() const {
        double den = RealPart * RealPart + ImaginaryPart * ImaginaryPart;
        if (den == 0) throw std::runtime_error("Division by zero");
        return {RealPart / den, -ImaginaryPart / den};
    }
    complex operator/(const complex& other) const { return (*this) * other.Invert(); }
    complex operator/(double val) const { 
        if (val == 0) throw std::runtime_error("Division by zero");
        return {RealPart / val, ImaginaryPart / val}; 
    }
    friend complex operator/(double val, const complex& c) { return complex(val) / c; }
//absolute value
    double absolute() const {
        return std::sqrt(RealPart * RealPart + ImaginaryPart * ImaginaryPart);
    }
//double conveersion
    operator double() const { return RealPart; }
//printing
    friend std::ostream& operator<<(std::ostream& os, const complex& c) {
        if (c.ImaginaryPart == 0) return os << c.RealPart;
        os << c.RealPart << (c.ImaginaryPart >= 0 ? " + " : " - ") 
           << std::abs(c.ImaginaryPart) << "i";
        return os;
    }
};

//Matrix struct - can handle: real and complex matrices - addition, subtraction, scalar mult 
//(compl, real), multiplication by matrix, vector multiplication, transpose, trace, 
template <typename T>
struct matrix {int rows, cols; std::vector<std::vector<T>> data;
    matrix(int r, int c, T val = T()) : rows(r), cols(c), data(r, std::vector<T>(c, val)) {}
    void fill(T value){
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                data[i][j] = value;
            }
            
        }
        
    }
    //KITERJESZTENI NEM NEGYZETESRE IS
    matrix<T> transpose(){

        matrix<T> result(cols, rows);

        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {   
               result.data[j][i] = data[i][j];
            }
        }
        return result;
    }
    matrix operator*(const matrix& other) const {
        matrix res(rows, other.cols, T()); // Initialize with "zero" of type T
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < other.cols; ++j)
                for (int k = 0; k < cols; ++k)
                    res.data[i][j] = res.data[i][j] + (data[i][k] * other.data[k][j]);
        return res;
    }
    matrix<T> operator*(const T& scalar) const {
        matrix<T> res(rows, cols);
        for (int i = 0; i < rows; ++i)
            for (int j = 0; j < cols; ++j)
                res.data[i][j] = data[i][j] * scalar;
        return res;
    }
    matrix<T> operator+(const matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for addition.");

        matrix<T> res(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                res.data[i][j] = data[i][j] + other.data[i][j];
            }
        }
        return res;
    }
    matrix<T> operator-(const matrix<T>& other) const {
        if (rows != other.rows || cols != other.cols)
            throw std::invalid_argument("Matrix dimensions must match for subtraction.");

        matrix<T> res(rows, cols);
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                res.data[i][j] = data[i][j] - other.data[i][j];
            }
        }
        return res;
    }
    friend matrix<T> operator*(const T& scalar, const matrix<T>& m) {
    return m * scalar; // Reuse the logic above
    } 
    std::vector<T> operator*(const std::vector<T>& vec) const {
        if (vec.size() != (size_t)cols)
            throw std::invalid_argument("Matrix columns must match vector size");

        std::vector<T> res(rows, T()); // Initialize result vector with zeros
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                res[i] = res[i] + (data[i][j] * vec[j]);
            }
        }
        return res;
}     
    matrix<complex> inverse() const{
        matrix<complex> result(2,2);
        if(rows==cols && rows==2){
            complex det = data[0][0]*data[1][1] - data[1][0]*data[0][1];
            if(det != 0){
            complex a = det.Invert();
            result.data[0][0] = a*data[1][1];
            result.data[0][1] = -a*data[0][1];
            result.data[1][0] = -a*data[1][0];
            result.data[1][1] = a*data[0][0];
            }
        }
        return result;
    }
    T trace(){
        T result = data[0][0];
        if (rows == cols){
            for (int i = 1; i < rows; i++)
            {
                result = result + data[i][i];
            }
        }
        else {
            std::cout << "Non-square matrix has no definiton for trace";
        }
        return result;
    }
    friend std::ostream& operator<<(std::ostream& os, const matrix<T>& m) {
        os << "\n";
        for (int i = 0; i < m.rows; ++i) {
            os << "|";
            for (int j = 0; j < m.cols; ++j) {
                os << m.data[i][j] << "\t"; // Uses T's << operator (e.g., Complex's)
            }
            os << "|";
            os << "\n";
            
        }
        return os;
    }

};

//temporary state struct, might not need it later, pretty unnecesary
//*********************************
struct state {matrix<complex> down; matrix<complex> up; complex alpha; complex beta;
    state(matrix<complex> d, matrix<complex> u, complex a = {1,0}, complex b = {1,0}): up(u), down(d), alpha(a), beta(b) {}

    state operator+(const state& other) const{
        state result({2,1}, {2,1});
        result.up.data = {{1,0},{0,0}};
        result.down.data = {{0,0}, {0,1}};
        result.alpha = alpha + other.alpha;
        result.beta = beta + other.beta;
        return result; 
    }
    state operator*(complex val){
        state result({2,1}, {2,1});
        result.up = val*alpha*up;
        result.down = val*beta*down;
        return result;
    }
    friend state operator*(const matrix<complex>& m, const state& s)
    {
        matrix<complex> placeholder(2,1);
        state result({2,1},{2,1});
        result.up.data = {{1,0}, {0,0}};
        result.down.data = {{0,0}, {0,1}};
        placeholder = s.alpha*(m*s.up) + s.beta*(m*s.down);
        result.alpha = placeholder.data[0][0];
        result.beta = placeholder.data[1][0];
        return result;
    }
    friend std::ostream& operator<<(std::ostream& os, const state& s){
        os << "( " << s.alpha << " ) |0> + ( " << s.beta << " ) |1>";
        return os;
    }

};
//*********************************


//Variables of projection
const int WIDTH = 240;
const int HEIGHT = 120;
const float PI = 3.1415926535f;

//Sphere properties
const float RADIUS = 20.0f;
const float ROTATION_SPEED = 0.00f;

//Helper vector struct for camera math
struct Vec3 {
    float x, y, z;
    Vec3 operator-(const Vec3& v) const { return {x - v.x, y - v.y, z - v.z}; }
    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }
    Vec3 cross(const Vec3& v) const { return {y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x}; }
    float norm() const { return std::sqrt(x * x + y * y + z * z); }
    Vec3 unit() const { float n = norm(); return {x / n, y / n, z / n}; }
};

//Draw function for the sphere and the axes, and the vector on the bloch sphere. indata is the
// vectors in  (theta, phi) coordinates
void Draw(std::vector<std::vector<double>> indata) {
    float angle = 0.0f;
    Vec3 cameraPos = {35.0f, 25.0f, 35.0f}; 
    Vec3 target = {0.0f, 0.0f, 0.0f}; 

    Vec3 forward = (target - cameraPos).unit();
    Vec3 tempUp = {0, 1, 0};
    Vec3 right = tempUp.cross(forward).unit();
    Vec3 up = forward.cross(right);

    float K1 = WIDTH * 0.8f;

    for (int k = 0; k < indata.size(); k++) {
        char buffer[HEIGHT][WIDTH];
        for (int i = 0; i < HEIGHT; ++i)
            for (int j = 0; j < WIDTH; ++j) buffer[i][j] = ' ';

        // --- 1. Draw Coordinate Axes (Fixed in World Space) ---
        float axisLen = RADIUS * 1.5f;
        Vec3 axes[3] = {{axisLen, 0, 0}, {0, axisLen, 0}, {0, 0, axisLen}};
        char axisLabels[3] = {'X', 'Y', 'Z'};

        for (int i = 0; i < 3; ++i) {
            for (float t = 0; t <= 1.0f; t += 0.02f) {
                Vec3 p = {axes[i].x * t, axes[i].y * t, axes[i].z * t};
                Vec3 rel = p - cameraPos;
                float zCam = rel.dot(forward);
                if (zCam > 0.1f) {
                    int xp = (int)(WIDTH / 2 + (K1 * rel.dot(right) / zCam));
                    int yp = (int)(HEIGHT / 2 - (K1 * 0.5f * rel.dot(up) / zCam));
                    if (xp >= 0 && xp < WIDTH && yp >= 0 && yp < HEIGHT) 
                        buffer[yp][xp] = (t > 0.95f) ? axisLabels[i] : '+';
                }
            }
        }

        // --- 2. Draw Sphere Mesh (Rotating) ---
        for (float v = 0; v <= PI; v += 0.05f) {
            for (float u = 0; u <= 2 * PI; u += 0.05f) {
                bool isEquator = std::abs(v - PI / 2.0f) < 0.025f;
                bool isLon1 = std::abs(u - 0.0f) < 0.025f || std::abs(u - PI) < 0.025f;
                bool isLon2 = std::abs(u - PI / 2.0f) < 0.025f || std::abs(u - 3 * PI / 2.0f) < 0.025f;

                if (isEquator || isLon1 || isLon2) {
                    float x = RADIUS * std::sin(v) * std::cos(u);
                    float y = RADIUS * std::sin(v) * std::sin(u);
                    float z = RADIUS * std::cos(v);

                    float rx = x * std::cos(angle) + z * std::sin(angle);
                    float rz = -x * std::sin(angle) + z * std::cos(angle);
                    
                    Vec3 rel = Vec3{rx, y, rz} - cameraPos;
                    float zCam = rel.dot(forward);
                    if (zCam > 0.1f) {
                        int xp = (int)(WIDTH / 2 + (K1 * rel.dot(right) / zCam));
                        int yp = (int)(HEIGHT / 2 - (K1 * 0.5f * rel.dot(up) / zCam));
                        if (xp >= 0 && xp < WIDTH && yp >= 0 && yp < HEIGHT) buffer[yp][xp] = '.';
                    }
                }
            }
        }

        // --- 3. Draw Silhouette (Circle orthogonal to projection plane) ---
        float distToCenter = cameraPos.norm();
        float screenRadius = (RADIUS * K1) / distToCenter; // Perspective scale
        for (float a = 0; a < 2 * PI; a += 0.02f) {
            int xp = (int)(WIDTH / 2 + screenRadius * std::cos(a));
            int yp = (int)(HEIGHT / 2 + screenRadius * 0.5f * std::sin(a)); // Aspect ratio fix
            if (xp >= 0 && xp < WIDTH && yp >= 0 && yp < HEIGHT) buffer[yp][xp] = '#';
        }

                // --- 4. Draw Specific Data Point & Vector from Origin ---
        float u_pt = (float)indata[k][0]; // theta
        float v_pt = (float)indata[k][1]; // phi

        // 3D Cartesian coordinates of the point
        float px = RADIUS * std::sin(v_pt) * std::cos(u_pt);
        float py = RADIUS * std::sin(v_pt) * std::sin(u_pt);
        float pz = RADIUS * std::cos(v_pt);

        // Rotate the point (sync with sphere rotation)
        float prx = px * std::cos(angle) + pz * std::sin(angle);
        float prz = -px * std::sin(angle) + pz * std::cos(angle);
        float pry = py;

        // Draw the line from (0,0,0) to the point
        for (float t = 0; t <= 1.0f; t += 0.05f) {
            Vec3 linePt = {prx * t, pry * t, prz * t};
            Vec3 lRel = linePt - cameraPos;
            float lzCam = lRel.dot(forward);
            
            if (lzCam > 0.1f) {
                int lx = (int)(WIDTH / 2 + (K1 * lRel.dot(right) / lzCam));
                int ly = (int)(HEIGHT / 2 - (K1 * 0.5f * lRel.dot(up) / lzCam));
                
                if (lx >= 0 && lx < WIDTH && ly >= 0 && ly < HEIGHT) {
                    // Use a different character for the vector line
                    if (buffer[ly][lx] == ' ') buffer[ly][lx] = '*'; 
                }
            }
        }

        // Draw the actual point at the end of the line
        Vec3 pRel = Vec3{prx, pry, prz} - cameraPos;
        float pzCam = pRel.dot(forward);
        if (pzCam > 0.1f) {
            int pxp = (int)(WIDTH / 2 + (K1 * pRel.dot(right) / pzCam));
            int pyp = (int)(HEIGHT / 2 - (K1 * 0.5f * pRel.dot(up) / pzCam));
            if (pxp >= 0 && pxp < WIDTH && pyp >= 0 && pyp < HEIGHT) {
                buffer[pyp][pxp] = 'O'; // The point itself
            }
        }


        // Render to console
        std::cout << "\x1b[H";
        for (int i = 0; i < HEIGHT; ++i) {
            for (int j = 0; j < WIDTH; ++j) std::cout << buffer[i][j];
            std::cout << "\n";
        }
        angle += ROTATION_SPEED;
        std::this_thread::sleep_for(std::chrono::milliseconds(80));
    }
}

//conversion from state to theta, phi - pretty unnecesary yet
std::vector<double> convert(state input)
{
    return {2*std::acos(input.alpha.absolute()),std::atan2(input.beta.ImaginaryPart, input.beta.RealPart) - std::atan2(input.alpha.ImaginaryPart, input.alpha.RealPart)};
}

int main(int argc, char* argv[])
{   

    //constants and shorthands
    complex id = {0,1};
    double hbar = 1;
    matrix<complex> cidentity(2,2);
    cidentity.data = {{{1,0},{0,0}},{{0,0},{1,0}}};

    //Pauli matrices
    matrix<complex> sx(2,2);
    matrix<complex> sy(2,2);
    matrix<complex> sz(2,2);

    sx.data = {{{0,0},{1,0}},{{1,0},{0,0}}};
    sy.data = {{{0,0},{0,-1}},{{0,1},{0,0}}};
    sz.data = {{{1,0},{0,0}},{{0,0},{-1,0}}};
    std::cout << sx << std::endl;
    std::cout << sy << std::endl;
    std::cout << sz << std::endl;

    std::cout << sx.inverse() <<std::endl;

    //Basis staetvectors
    matrix<complex>  up0(2,1);
    up0.data = {{1,0}, {0,0}};
    matrix<complex>  down0(2,1);
    down0.data = {{0,0}, {1,0}};

    //state??
    state state0({2,1}, {2,1});
    state0.up = up0; 
    state0.down = down0;
    state0.alpha = (1/sqrt(2));
    state0.beta = (1/sqrt(2));



    /*
    std::cout << sx.data[0][1];

    
    matrix<complex> H = omegax*sx + omegay*sy + omegaz*sz;
    double deltat = 0.001;
    matrix<complex> U(2,2);
    std::cout << H << std::endl;
    std::cout << state0 << std::endl;
    std::vector<std::vector<double>> data0;
    for (int i = 0; i < 30; i++)
    {
        U = cidentity - ((i/hbar)*deltat)*H;
        state0 = U*state0;
        std::vector<double> converted = convert(state0);
        std::cout << "Thtea = " << converted[0] << " Phi = " << converted[1] << std::endl;
        data0.push_back({converted[0], converted[1]});

    }
    
    matrix<double> to_transpose(3,2);
    to_transpose.data = {{{1,2}},{{4,5}},{{7,8}}};
    std::cout << to_transpose << std::endl;
    std::cout << to_transpose.transpose() << std::endl;
    Draw(data0);
    */
    
    //Definition of the Hamiltonian
    complex omegax = {3,1};
    complex omegay = {3,2};
    complex omegaz = {5,8};
    std::cout << state0 << std::endl;
    matrix<complex> H(2,2);
    H = ((hbar*omegaz)/2.0)*sz;

    //the statevector from the state struct - still redundant
    matrix<complex> statevector = state0.alpha*state0.up + state0.beta*state0.down;

    //density matrix
    matrix<complex> rho = statevector*(statevector.transpose());
    matrix<complex> rhot(2,2);
    double t = 0.0;
    complex rhot12;
    complex rhot21;
    std::vector<std::vector<double>> save;
    
    for (int i = 0; i < 130; i++)
    {
        rhot12 = {0.5*std::cos(-omegaz*t), 0.5*std::sin(-omegaz*t)};
        rhot21 = {0.5*std::cos(omegaz*t), 0.5*std::sin(omegaz*t)};
        rhot.data = {{{0.5, 0},rhot12},{rhot12,{0.5}}};


        complex x = rhot12 + rhot21;
        complex y = id*(rhot12-rhot21);
        complex z = rhot.data[0][0] - rhot.data[1][1];

        double radius = std::sqrt(double(x)*double(x) + double(y)*double(y) + double(z)*double(z));
        double theta = std::atan2(double(y), double(x));
        double phi = std::acos(double(z)/radius);
        save.push_back({theta, phi, radius});

        t += 0.01;
    }
    Draw(save);
    
}