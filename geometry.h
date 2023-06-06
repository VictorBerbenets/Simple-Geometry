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
                                                //  ------------------
class Line {                                    //  |ax + by + c = 0 |  
private:                                        //  ------------------                         
    std::vector<double> dir_vect_;
    double koeff_  = 0; // = k
    double offset_ = 0; // = b
public:
    Line(double koeff, double offset): koeff_(koeff), offset_(offset) {};
    Line(const Point& point1, const Point& point2);
    Line(const Point& point, double koeff);
    Line(double a, double b, double c);
    Line(const Line& line) = default;

    bool operator==(const Line& line) const;
    bool operator!=(const Line& line) const;

    double Getk() const {
        return koeff_;
    }

    double Getb() const {
        return offset_;
    }

    double Getv(size_t id) const {
        return dir_vect_[id];
    }

    ~Line() = default;
};

Line::Line(double a, double b, double c) {
    dir_vect_.reserve(3);

    dir_vect_.push_back(-b);
    dir_vect_.push_back(a);
    dir_vect_.push_back(c);
}

Line::Line(const Point& point1, const Point& point2) {
    double x_offset = point2.Getx() - point1.Getx();
    double y_offset = point2.Gety() - point1.Gety();

    dir_vect_.reserve(2);
    dir_vect_.push_back(x_offset);  // line vector
    dir_vect_.push_back(y_offset);

    if (!x_offset) {
        koeff_ = std::numeric_limits<double>::infinity();
    } else {
        koeff_ = y_offset / x_offset;
    }

    offset_ = point1.Gety() - koeff_ * point1.Getx();
}

Line::Line(const Point& point, double koeff): koeff_(koeff), offset_(point.Gety() - koeff_ * point.Getx()) {
    dir_vect_.push_back(1);
    dir_vect_.push_back(koeff_);
}

bool Line::operator==(const Line& line) const {
    return is_equal::is_equal(koeff_, line.koeff_) && is_equal::is_equal(offset_, line.offset_);
}

bool Line::operator!=(const Line& line) const {
    return !(*this == line);
}


class Shape {
    virtual double area()     const = 0;
    virtual double perimetr() const = 0;
    
    virtual bool operator==(const Shape& another)    const = 0;
    virtual bool isCongruentTo(const Shape& another) const = 0;
    virtual bool isSimilarTo(const Shape& another)   const = 0;
    virtual bool containsPoint(const Point& another) const = 0;

    virtual void rotate(const Point& center, double angle)      = 0;
    virtual void scale(const Point& center, double coefficient) = 0;
    virtual void reflect(const Point& center) = 0;
    virtual void reflect(const Line& axis)    = 0;
};

//=============================================Class Poligon======================================================//
class Poligon: public Shape {
private:
    double mul_vect(const Line& line1, const Line& line2) const; //psevdo vector product
private:
    std::vector<Point> vertices_;
public:
    /// 
    Poligon(const std::vector<Point>& vertices);
    // Poligon(const Point&...);
    Poligon(const Poligon& copy) = default;
    ///  

    size_t verticesCount() const;
    std::vector<Point> getVertices() const;
    bool isConvex() const;

//virtual funcs
    double area()     const override;
    double perimetr() const override;
    
    bool operator==(const Shape& another)    const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another)   const override;
    bool containsPoint(const Point& point)   const override;

    void rotate(const Point& center, double angle)      override;
    void scale(const Point& center, double coefficient) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis)    override;
//
    ~Poligon() = default;
};

Poligon::Poligon(const std::vector<Point>& vertices) {
    size_t copy_size = vertices.size();
    vertices_.reserve(copy_size);

    for (size_t current_id = 0; current_id < copy_size; ++current_id) {
        vertices_.push_back(vertices[current_id]);
    }
}


size_t Poligon::verticesCount() const {
    return vertices_.size();
}

std::vector<Point> Poligon::getVertices() const {
    return vertices_;
}

double Poligon::mul_vect(const Line& line1, const Line& line2) const {                    //   |a1 a2| 
    return line1.Getv(0) * line2.Getv(1) - line1.Getv(1) * line2.Getv(0);                 //   |b1 b2|  = a1*b2 - a2*b1
}                                                       

bool Poligon::isConvex() const {
    size_t limit = vertices_.size();

    for (size_t id = 0; id < limit; ++id) {
        Line line1(vertices_.at(id), vertices_.at( (id + 1) % limit) );
        Line line2(vertices_.at( (id + 1) % limit ), vertices_.at( (id + 2) % limit) );
        Line line3(vertices_.at( (id + 2) % limit ), vertices_.at( (id + 3) % limit) );

        // std::cout << "line1 * line2 = " << mul_vect(line1, line2) << std::endl;
        // std::cout << "line2 * line3 = " << mul_vect(line2, line3) << std::endl;
        // std::cout << "line1: " << " (" << line1.Getv(0) << " ; " << line1.Getv(1) << ")" << std::endl;
        // std::cout << "line2: " << " (" << line2.Getv(0) << " ; " << line2.Getv(1) << ")" << std::endl;
        // std::cout << "line3: " << " (" << line3.Getv(0) << " ; " << line3.Getv(1) << ")" << std::endl;

        if (mul_vect(line1, line2) * mul_vect(line2, line3) < 0) {
            return false;
        }
        std::cout << '\n';
    }
    return true;
}


double Poligon::area() const {

}

double Poligon::perimetr() const {
    
}

bool Poligon::operator==(const Shape& another) const {

}

bool Poligon::isCongruentTo(const Shape& another) const {

}

bool Poligon::isSimilarTo(const Shape& another) const {

}

bool Poligon::containsPoint(const Point& point) const {

}

void Poligon::rotate(const Point& center, double angle) {

}

void Poligon::scale(const Point& center, double angle) {
    
}

void Poligon::reflect(const Point& center) {
    
}

void Poligon::reflect(const Line& axis) {
    
}
//================================================================================================================//


//================================================================================================================//
class Ellipse: public Shape {
private:
    Point focus1_;
    Point focus2_;
    double a_ = 0;
    double b_ = 0;
    double e_ = 0;    //eccentricity
public:

    Ellipse(const Point& point1, const Point& point2, double dist);

    std::pair<Point, Point> focuses() const;
    std::pair<Line, Line> directrices() const;
    double eccentricity() const;
    Point center() const;

//virtual funcs
    double area()     const override;
    double perimetr() const override;
    
    bool operator==(const Shape& another)    const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another)   const override;
    bool containsPoint(const Point& point)   const override;

    void rotate(const Point& center, double angle)      override;
    void scale(const Point& center, double coefficient) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis)    override;
//

    virtual ~Ellipse() = default;
};

Ellipse::Ellipse(const Point& point1, const Point& point2, double dist): focus1_(point1), focus2_(point2) {
    a_ = dist / 2;
    b_ = sqrt( pow(focus1_.Getx() - focus2_.Getx(), 2) + pow(focus1_.Gety() - focus2_.Gety(), 2) );
    e_ = sqrt( pow(a_, 2) - pow(b_, 2));
}

double Ellipse::eccentricity() const {
    return e_;
}

Point Ellipse::center() const {
    Point ret_value( (focus1_.Getx() + focus2_.Getx()) / 2, (focus1_.Gety() + focus2_.Gety()) / 2);
    return ret_value;
}


std::pair<Line, Line> Ellipse::directrices() const {
    Line direct_1(0, 1, -(a_ / e_));
    Line direct_2(0, 1, a_ / e_);

    std::pair<Line, Line> ret_value(direct_1, direct_2);
    return ret_value;
}

std::pair<Point, Point> Ellipse::focuses() const {
    std::pair<Point, Point> focuses(focus1_, focus2_);
    return focuses;
}




double Ellipse::area() const {

}

double Ellipse::perimetr() const {
    
}

bool Ellipse::operator==(const Shape& another) const {

}

bool Ellipse::isCongruentTo(const Shape& another) const {

}

bool Ellipse::isSimilarTo(const Shape& another) const {

}

bool Ellipse::containsPoint(const Point& point) const {

}

void Ellipse::rotate(const Point& center, double angle) {

}

void Ellipse::scale(const Point& center, double angle) {
    
}

void Ellipse::reflect(const Point& center) {
    
}

void Ellipse::reflect(const Line& axis) {
    
}
//================================================================================================================//

//================================================================================================================//
class Circle: public Ellipse {
    using Ellipse::Ellipse;
private:
    double r_;
    Point center_;
public:
    Circle(const Point& center, double radius): Ellipse(center, center, radius), center_(center), r_(radius) {}
//virtual funcs
    double area()     const override;
    double perimetr() const override;
    
    bool operator==(const Shape& another)    const override;
    bool isCongruentTo(const Shape& another) const override;
    bool isSimilarTo(const Shape& another)   const override;
    bool containsPoint(const Point& point)   const override;

    void rotate(const Point& center, double angle)      override;
    void scale(const Point& center, double coefficient) override;
    void reflect(const Point& center) override;
    void reflect(const Line& axis)    override;
//
    ~Circle() = default;
};


double Circle::area() const {
    return M_PI * pow(r_, 2);
}

double Circle::perimetr() const {
    return 2 * M_PI * r_;
}

bool Circle::operator==(const Shape& another) const {

}

bool Circle::isCongruentTo(const Shape& another) const {

}

bool Circle::isSimilarTo(const Shape& another) const {

}

bool Circle::containsPoint(const Point& point) const {

}

void Circle::rotate(const Point& center, double angle) {

}

void Circle::scale(const Point& center, double angle) {
    
}

void Circle::reflect(const Point& center) {
    
}

void Circle::reflect(const Line& axis) {
    
}
//================================================================================================================//

//================================================================================================================//
class Rectangle: public Poligon {
    using Poligon::Poligon;
private:
    Point left_;
    Point right_;

public:
    // Rectangle(const Point& left, const Point& right):  left_(left), right_(right){}
    Point center() const;
    std::pair<Line, Line> diagonals() const;



};

//================================================================================================================//

#endif