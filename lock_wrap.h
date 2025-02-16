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

	T *lock() const {
		void *p = (void *)&mutex_;
		((MUTEX_T *)p)->lock();
		return obj;
	}

	void unlock() const {
		void *p = (void *)&mutex_;
		((MUTEX_T *)p)->unlock();
	}
};
#endif