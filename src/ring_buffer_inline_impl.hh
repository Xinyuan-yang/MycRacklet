/* -------------------------------------------------------------------------- */
template<typename T>
inline RingBuffer<T>::RingBuffer(){

  values = NULL;
  size = 0;
  n = 0;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline RingBuffer<T>::~RingBuffer(){

}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::init(T * field_val, UInt size, T default_val){

  this->size = size;
  this->values = field_val;
  n=0;
  T * it = values;
  for (UInt i = 0; i < size; ++i) {
    *it = default_val;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::operator<<(T new_val){

  *(values+n%size) = new_val;
  ++n;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::current() const {
  return values+(n-1)%size;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::begin() const {
  return values;
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline T * RingBuffer<T>::end() const {
  return values+std::min(size-1,n-1);
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline void RingBuffer<T>::printself(std::ostream & stream) const {

  T * it = values;
  for (UInt i = 0; i < size; ++i) {
    stream << " " << *it;
    ++it;
  }
}

/* -------------------------------------------------------------------------- */
template<typename T>
inline std::ostream & operator <<(std::ostream & stream, const RingBuffer<T> & _this)
{
  _this.printself(stream);
  return stream;
}
