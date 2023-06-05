#ifndef geometry_h
#define geometry_h

#include <iostream>
#include <vector>
#include <cmath>

namespace is_equal{
    const double Epsilon = 1e-17;
    bool is_equal(double first_value, double second_value);
}

bool is_equal::is_equal(double first_value, double second_value) {
    return fabs(first_value - second_value) < Epsilon;
}

class Point {
private:
    double x_;
    double y_;
public:
    Point() = default;
    Point(double x, double y): x_(x), y_(y) {};
    Point(const Point&) = default;

    bool operator==(const Point& another) const{
        return is_equal::is_equal(x_, another.x_) && is_equal::is_equal(y_, another.y_);
    }

    bool operator!=(const Point& another) const{
        return !(*this == another);
    }

    double Getx() const {
        return x_;
    }
    double Gety() const {
        return y_;
    }

    ~Point() = default;
};

class Line {
private:
    double koeff_  = 0;
    double offset_ = 0;
public:
    Line(double koeff, double offset): koeff_(koeff), offset_(offset) {};
    Line(const Point& point1, const Point& point2);
    Line(const Point& point, double koeff);
    Line(const Line& line) = default;

    bool operator==(const Line& line);
    bool operator!=(const Line& line);

    double Get_Koeff() {
        return koeff_;
    }
    double Get_Offset() {
        return offset_;
    }

    ~Line() = default;
};

Line::Line(const Point& point1, const Point& point2) {
    double x_offset = point2.Getx() - point1.Getx();
    double y_offset = point2.Gety() - point2.Gety();

    if (!x_offset) {
        koeff_ = 0;
    } else {
        koeff_ = y_offset / x_offset;
    }

    offset_ = koeff_ * point1.Getx() - point1.Gety();
}

Line::Line(const Point& point, double koeff): koeff_(koeff), offset_(koeff_ * point.Getx() - point.Gety()) {}

bool Line::operator==(const Line& line) {
    return is_equal::is_equal(koeff_, line.koeff_) && is_equal::is_equal(offset_, line.offset_);
}

bool Line::operator!=(const Line& line) {
    return !(*this == line);
}


class Shape {
    virtual double area()     = 0;
    virtual double perimetr() = 0;
    
    virtual bool operator==(const Shape& another)    = 0;
    virtual bool isCongruentTo(const Shape& another) = 0;
    virtual bool isSimilarTo(const Shape& another)   = 0;
    virtual bool containsPoint(const Point& another) = 0;

    virtual void rotate(const Point& center, double angle)      = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis)    = 0;
};


class Poligon: public Shape {
private:
    std::vector<Point> vertices_;
public:
    Poligon(const std::vector<Point>& vertices);
    // Poligon(const Point&...);
    Poligon(const Poligon& copy) = default;
    
    size_t verticesCount();
    std::vector<Point> getVertices();
    bool isConvex();


    ~Poligon() = default;
};

class Ellipse: public Shape {
private:
    Point focus1_;
    Point focus2_;
    double a_ = 0;
    double b_ = 0;
    double eccentricity_;
public:

    Ellipse(const Point& point1, const Point& point2, double dist);


    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

    ~Ellipse() = default;
};

Ellipse::Ellipse(const Point& point1, const Point& point2, double dist): focus1_(point1), focus2_(point2) {
    a_ = dist / 2;
    b_ = sqrt( pow(focus1_.Getx() - focus2_.Getx(), 2) + pow(focus1_.Gety() - focus2_.Gety(), 2) );
    eccentricity_ = sqrt( pow(a_, 2) - pow(b_, 2));
}

double Ellipse::eccentricity() const {
    return eccentricity_;
}

Point Ellipse::center() const {
    Point ret_value( (focus1_.Getx() + focus2_.Getx()) / 2, (focus1_.Gety() + focus2_.Gety()) / 2);
    return ret_value;
}


#endif