#ifndef LOCK_WRAP_H
#define LOCK_WRAP_H

#include "threading.h"

template <typename T>
class LockWrap {
  private:
	T *obj;
	MUTEX_T mutex_;

  public:
	LockWrap(T *obj) : obj(obj) {};
	~LockWrap() { delete obj; };

	T *lock() {
		mutex_.lock();
		return obj;
	}

	void unlock() {
		mutex_.unlock();
	}
};
#endif