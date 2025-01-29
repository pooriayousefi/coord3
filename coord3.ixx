export module coord3;

import <stdexcept>;
import <type_traits>;
import <concepts>;
import <limits>;
import <cmath>;
import <tuple>;
import <iostream>;
import <fstream>;

// arithmetic concept definition
export template<class T> concept arithmetic = std::integral<T> || std::floating_point<T>;

// 3D space
export namespace R3
{
	// Cartesian coordinate system
	template<arithmetic T> struct cartesian_coordinate
	{
		// type definitions
		using value_type = T;
		using type = cartesian_coordinate<value_type>;

		// aligned constant member
		alignas(value_type) const value_type x, y, z;

		// default constructor (deleted)
		cartesian_coordinate() = delete;

		// constructor
		cartesian_coordinate(value_type X, value_type Y, value_type Z) :
			x{ static_cast<const value_type>(X) }, 
			y{ static_cast<const value_type>(Y) },
			z{ static_cast<const value_type>(Z) }
		{
		}
	};

	// arithmetic operator overloadings
	template<arithmetic T> constexpr bool operator==(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		if constexpr (std::is_floating_point_v<T>)
			return std::abs(lhs.x - rhs.x) < std::numeric_limits<T>::epsilon() &&
			std::abs(lhs.y - rhs.y) < std::numeric_limits<T>::epsilon() &&
			std::abs(lhs.z - rhs.z) < std::numeric_limits<T>::epsilon();
		else if constexpr (std::is_integral_v<T>)
			return lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z;
	}
	template<arithmetic T> constexpr bool operator!=(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		if constexpr (std::is_floating_point_v<T>)
			return std::abs(lhs.x - rhs.x) >= std::numeric_limits<T>::epsilon() &&
			std::abs(lhs.y - rhs.y) >= std::numeric_limits<T>::epsilon() &&
			std::abs(lhs.z - rhs.z) >= std::numeric_limits<T>::epsilon();
		else if constexpr (std::is_integral_v<T>)
			return lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z;
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator+(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		return { lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z };
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator-(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		return { lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z };
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator*(const cartesian_coordinate<T>& lhs, T scalar)
	{
		return { lhs.x * scalar, lhs.y * scalar, lhs.z * scalar };
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator*(T scalar, const cartesian_coordinate<T>& rhs)
	{
		return rhs * scalar;
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator/(const cartesian_coordinate<T>& lhs, T scalar)
	{
		if constexpr (std::is_floating_point_v<T>)
		{
			if (std::abs(scalar) < std::numeric_limits<T>::epsilon)
				throw std::overflow_error("Division by zero");
			return { lhs.x / scalar, lhs.y / scalar, lhs.z / scalar };
		}
		else if constexpr (std::is_integral_v<T>)
		{
			if (scalar == 0)
				throw std::overflow_error("Division by zero");
			return { lhs.x / scalar, lhs.y / scalar, lhs.z / scalar };
		}
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator/(T scalar, const cartesian_coordinate<T>& rhs)
	{
		if constexpr (std::is_floating_point_v<T>)
		{
			if (
				rhs.x < std::numeric_limits<T>::epsilon || 
				rhs.y < std::numeric_limits<T>::epsilon ||
				rhs.z < std::numeric_limits<T>::epsilon
				)
				throw std::overflow_error("Division by zero");
			return { scalar / rhs.x, scalar / rhs.y, scalar / rhs.z };
		}
		else if constexpr (std::is_integral_v<T>)
		{
			if (rhs.x == 0 || rhs.y == 0 || rhs.z == 0)
				throw std::overflow_error("Division by zero");
			return { scalar / rhs.x, scalar / rhs.y, scalar / rhs.z };
		}
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> operator-(const cartesian_coordinate<T>& rhs)
	{
		return { -rhs.x, -rhs.y, -rhs.z };
	}
	
	// apply meta-function
	template<arithmetic T, class F>
	constexpr cartesian_coordinate<T> apply(F&& f, const cartesian_coordinate<T>& c)
	{
		if constexpr (std::is_invocable_r_v<cartesian_coordinate<T>, F, const cartesian_coordinate<T>&>)
		{
			return f(c);
		}
		else if constexpr (std::is_invocable_r_v<T, F, const T&>)
		{
			return { f(c.x), f(c.y), f(c.z) };
		}
		else
		{
			static_assert(false, "Invalid function type");
		}
	}

	// vectorial operator overloadings
	template<arithmetic T> constexpr T dot(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> cross(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
	{
		return { lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x };
	}
	
	// metrics
	// Minkowski metric
	namespace Minkowski
	{
		template<arithmetic T> constexpr T distance(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs, T p)
		{
			return std::pow(std::pow(std::abs(lhs.x - rhs.x), p) + std::pow(std::abs(lhs.y - rhs.y), p) + std::pow(std::abs(lhs.z - rhs.z), p), 1 / p);
		}
		template<arithmetic T> constexpr T norm(const cartesian_coordinate<T>& rhs, T p)
		{
			return std::pow(std::pow(std::abs(rhs.x), p) + std::pow(std::abs(rhs.y), p) + std::pow(std::abs(rhs.z), p), 1 / p);
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> normalize(const cartesian_coordinate<T>& rhs, T p)
		{
			const auto n = norm(rhs, p);
			if (n < std::numeric_limits<T>::epsilon)
				throw std::overflow_error("Division by zero");
			return rhs / n;
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> direction(const cartesian_coordinate<T>& from, const cartesian_coordinate<T>& to, T p)
		{
			return normalize(to - from, p);
		}
	}

	// Manhattan metric (L1 norm)
	namespace Manhattan
	{
		template<arithmetic T> constexpr T distance(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
		{
			return std::abs(lhs.x - rhs.x) + std::abs(lhs.y - rhs.y) + std::abs(lhs.z - rhs.z);
		}
		template<arithmetic T> constexpr T norm(const cartesian_coordinate<T>& rhs)
		{
			return std::abs(rhs.x) + std::abs(rhs.y) + std::abs(rhs.z);
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> normalize(const cartesian_coordinate<T>& rhs)
		{
			const auto n = norm(rhs);
			if (n < std::numeric_limits<T>::epsilon)
				throw std::overflow_error("Division by zero");
			return rhs / n;
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> direction(const cartesian_coordinate<T>& from, const cartesian_coordinate<T>& to)
		{
			return normalize(to - from);
		}
	}

	// Euclidean metric (L2 norm)
	namespace Euclidean
	{
		template<arithmetic T> constexpr T distance(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
		{
			return std::sqrt(
				((lhs.x - rhs.x)* (lhs.x - rhs.x)) + 
				((lhs.y - rhs.y)* (lhs.y - rhs.y)) +
				((lhs.z - rhs.z)* (lhs.z - rhs.z))
			);
		}
		template<arithmetic T> constexpr T norm(const cartesian_coordinate<T>& rhs)
		{
			return std::sqrt(rhs.x * rhs.x + rhs.y * rhs.y + rhs.z * rhs.z);
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> normalize(const cartesian_coordinate<T>& rhs)
		{
			const auto n = norm(rhs);
			if (n < std::numeric_limits<T>::epsilon)
				throw std::overflow_error("Division by zero");
			return rhs / n;
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> direction(const cartesian_coordinate<T>& from, const cartesian_coordinate<T>& to)
		{
			return normalize(to - from);
		}
		template<arithmetic T> constexpr T Cosine_similarity(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
		{
			return dot(lhs, rhs) / (norm(lhs) * norm(rhs));
		}
	}

	// Chebyshev metric (Linf norm)
	namespace Chebyshev
	{
		template<arithmetic T> constexpr T distance(const cartesian_coordinate<T>& lhs, const cartesian_coordinate<T>& rhs)
		{
			return std::max({ std::abs(lhs.x - rhs.x), std::abs(lhs.y - rhs.y), std::abs(lhs.z - rhs.z) });
		}
		template<arithmetic T> constexpr T norm(const cartesian_coordinate<T>& rhs)
		{
			return std::max({ std::abs(rhs.x), std::abs(rhs.y), std::abs(rhs.z) });
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> normalize(const cartesian_coordinate<T>& rhs)
		{
			const auto n = norm(rhs);
			if (n < std::numeric_limits<T>::epsilon)
				throw std::overflow_error("Division by zero");
			return rhs / n;
		}
		template<arithmetic T> constexpr cartesian_coordinate<T> direction(const cartesian_coordinate<T>& from, const cartesian_coordinate<T>& to)
		{
			return normalize(to - from);
		}
	}

	// to curvilinear coordinate systems
	template<arithmetic T> constexpr std::tuple<T, T, T> to_cylindrical_coordinate(const cartesian_coordinate<T>& c)
	{
		return {
			std::hypot(c.x, c.y),
			std::atan2(c.y, c.x),
			c.z
		};
	}
	template<arithmetic T> constexpr std::tuple<T, T, T> to_spherical_coordinate(const cartesian_coordinate<T>& c)
	{
		return {
			Euclidean::norm(c),
			std::acos(c.z / Euclidean::norm(c)),
			std::atan2(c.y, c.x)
		};
	}

	// from curvilinear coordinate systems
	template<arithmetic T> constexpr cartesian_coordinate<T> from_cylindrical_coordinate(T rho, T phi, T Z)
	{
		return { rho * std::cos(phi), rho * std::sin(phi), Z };
	}
	template<arithmetic T> constexpr cartesian_coordinate<T> from_spherical_coordinate(T r, T theta, T phi)
	{
		return { r * std::sin(theta) * std::cos(phi), r * std::sin(theta) * std::sin(phi), r * std::cos(theta) };
	}

	// stream operator overloadings
	template<std::floating_point T> constexpr std::ostream& operator<<(std::ostream& os, const cartesian_coordinate<T>& c)
	{
		os << '(' << c.x << ',' << c.y << ',' << c.z << ')';
		return os;
	}
	template<std::floating_point T> constexpr std::ofstream& operator<<(std::ofstream& os, const cartesian_coordinate<T>& c)
	{
		os << '(' << c.x << ',' << c.y << ',' << c.z << ')';
		return os;
	}
}
