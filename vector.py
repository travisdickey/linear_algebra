## class containing various vector operations

from math import sqrt, acos, pi
from decimal import Decimal, getcontext

getcontext().prec = 10

class Vector(object):
    '''takes coordinates of a vector; contains various functions
       for vector operations'''

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG = 'No unique orthogonal component'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component'
    ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG = 'Only defined in two or three dimensions'

    def __init__(self, coordinates):
        '''set coordinates as tuple'''
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def is_orthogonal_to(self, v, tolerance=1e-10):
        '''takes two vectors and returns boolean on orthogonality'''
        return abs(self.dot(v)) < tolerance

    def is_parallel_to(self, v):
        '''takes two vectors and returns boolean on parallelism'''
        return ( self.is_zero() or
                v.is_zero() or
                self.angle_with(v) == 0 or
                self.angle_with(v) == pi )

    def is_zero(self, tolerance=1e-10):
        '''takes vector and checks whether magnitude is 0'''
        return self.magnitude() < tolerance

    def magnitude(self):
        '''takes vector and calculates magnitude'''
        squared_coordinates = [x**2 for x in self.coordinates]
        return sqrt(sum(squared_coordinates))

    def normalized(self):
        '''takes vector and returns its unit vector'''
        try:
            magnitude = self.magnitude()
            return Vector(self.times_scalar(Decimal('1.0') / Decimal(magnitude)))
        except ZeroDivisionError:
            raise Exception(CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot(self, v):
        '''takes two vectors and returns dot product'''
        return sum([x * y for x,y in zip(self.coordinates, v.coordinates)])

    def area_of_triangle_with(self, v):
        '''takes two vectors and returns area of triangle formed'''
        return Decimal(self.area_of_paralellogram_with(v)) / Decimal('2.0')

    def area_of_paralellogram_with(self, v):
        '''takes two vectors and returns area of parallelogram formed'''
        cross_product = self.cross(v)
        return cross_product.magnitude()

    def cross(self, v):
        '''takes two vectors and returns cross product'''
        try:
            x_1, y_1, z_1 = self.coordinates
            x_2, y_2, z_2 = v.coordinates
            new_coordinates = [ y_1*z_2 - y_2*z_1 ,
                                -(x_1*z_2 - x_2*z_1),
                                x_1*y_2 - x_2*y_1   ]
            return Vector(new_coordinates)

        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                self_embedded_in_R3 = Vector(self.coordinates + ('0',))
                v_embedded_in_R3 = vector(v.coordinates + ('0',))
                return self_embedded_in_R3.cross(v_embedded_in_R3)
            elif (msg == 'too many values to unpack' or
                  msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINED_IN_TWO_THREE_DIMS_MSG)
            else:
                raise e

    def component_orthogonal_to(self, basis):
        '''takes two vectors and returns the orthogonal component'''
        try:
            projection = Vector(self.component_parallel_to(basis))
            return self.minus(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_ORTHOGONAL_COMPONENT_MSG)
            else:
                raise e

    def component_parallel_to(self, basis):
        '''takes two vectors and returns the parallel component'''
        try:
            u = basis.normalized()
            weight = self.dot(u)
            return u.times_scalar(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def angle_with(self, v, in_degrees=False):
        '''takes two vectors and returns the angle formed by the two'''
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            angle_in_radians = acos(round(float(u1.dot(u2)),4))

            if in_degrees:
                degrees_per_radian = 180. / pi
                return angle_in_radians * degrees_per_radian
            else:
                return angle_in_radians
        except Exception as e:
            print e
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')

    def plus(self, v):
        '''takes two vectors and adds them'''
        new_coordinates = [x + y for x,y in zip(self.coordinates, v.coordinates)]
        return new_coordinates

    def minus(self, v):
        '''takes two vectors and substracts them'''
        new_coordinates = [x - y for x,y in zip(self.coordinates, v.coordinates)]
        return new_coordinates

    def times_scalar(self, c):
        '''takes a vector and scalar and multiplies them'''
        new_coordinates = [Decimal(c)*x for x in self.coordinates]
        return new_coordinates

    def __getitem__(self, i):
        return self.coordinates[i]

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        '''takes two vectors and returns boolean on equality'''
        return self.coordinates == v.coordinates
