#ifndef TYPE_NAME_HPP
#define TYPE_NAME_HPP


#include <string>


template <typename T>
class TypeName {
public:
	static inline std::string getName()
	{
		return typeid(T).name();
	}
};


template <>
class TypeName<long> {
public:
	static inline std::string getName()
	{
		return "long";
	}
};


template <>
class TypeName<unsigned long> {
public:
	static inline std::string getName()
	{
		return "unsigned long";
	}
};


template <>
class TypeName<int> {
public:
	static inline std::string getName()
	{
		return "int";
	}
};


template <>
class TypeName<unsigned int> {
public:
	static inline std::string getName()
	{
		return "unsigned int";
	}
};


template <>
class TypeName<short> {
public:
	static inline std::string getName()
	{
		return "short";
	}
};


template <>
class TypeName<unsigned short> {
public:
	static inline std::string getName()
	{
		return "unsigned short";
	}
};


template <>
class TypeName<char> {
public:
	static inline std::string getName()
	{
		return "char";
	}
};


template <>
class TypeName<unsigned char> {
public:
	static inline std::string getName()
	{
		return "unsigned char";
	}
};


template <>
class TypeName<std::string> {
public:
	static inline std::string getName()
	{
		return "string";
	}
};


template <>
class TypeName<bool> {
public:
	static inline std::string getName()
	{
		return "bool";
	}
};


template <>
class TypeName<float> {
public:
	static inline std::string getName()
	{
		return "float";
	}
};


template <>
class TypeName<double> {
public:
	static inline std::string getName()
	{
		return "double";
	}
};


#endif // TYPE_NAME_HPP
