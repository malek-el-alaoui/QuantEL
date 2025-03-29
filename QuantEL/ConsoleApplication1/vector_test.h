#pragma once
#include <iostream>
using namespace std;

template <typename T> class vector_test
{
  private:
	  T* arr;
	  int size;
	  int capacity;
  public:
	  vector_test() : arr(nullptr), size(0), capacity(0) {};
	  ~vector_test(){ delete[] arr }
	  void push_back(const T& x) 
	  {
		  if (size < capacity)
		  {
			  arr[size++] = x;
		  }
		  else
		  {
			  delete[] arr;
			  T* new_data = new T[2 * capacity];
			  for (int i = 0; i < size + 1; ++i) {
				  new_data = data[i];
			  }
			  new_data[size] = x
		  }
		  arr = new_data;
	  };

};

