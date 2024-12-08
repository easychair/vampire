#pragma once

#include <vector>

#include <iostream>

namespace Lib {

  template <typename T>
  class SparseMatrix {
    private:
      class Slice {
        friend class SparseMatrix<T>;
        public:
          Slice() = default;

          void set(unsigned index, T value);

          /**
           * @brief Returns a pointer to the value at the given index. If no value is set at the given index, returns nullptr.
           * @note This function is not marked const because it modifies the lastCheckedIndex and lastFoundIndex members.
           */
          T get(unsigned index, T def) const;

          void del(unsigned index);

          unsigned size() const;

        private:
          std::vector<std::pair<unsigned, T>> data;

          /**
           * @note This function is not marked const because it modifies the lastCheckedIndex and lastFoundIndex members.
           */
          std::pair<bool, unsigned> find(unsigned index) const;

          /**
           * @note Those fields are transparent to the user, so they are mutable.
           */
          mutable unsigned lastCheckedIndex = 0xFFFFFFFF;
          mutable unsigned lastFoundIndex = 0xFFFFFFFF;

      };
    public:
      SparseMatrix(unsigned rows, unsigned cols);

      void set(unsigned row, unsigned col, T value);

      void del(unsigned row, unsigned col);

      /**
       * @brief Returns a pointer to the value at the given index. If no value is set at the given index, returns nullptr.
       */
      T get(unsigned row, unsigned col, T def) const;

      /**
       * @brief Returns a reference to the non-zero elements in the given row.
       */
      std::vector<std::pair<unsigned, T>>& getSetOnRow(unsigned row);

      void swapRows(unsigned row1, unsigned row2);

      void reset();

      void reshape(unsigned rows, unsigned cols);

    private:
      std::vector<Slice> data;
      unsigned rows;
      unsigned cols;

  };

  /***********************************************************************************************/
  /*                                            SLICE                                            */
  /***********************************************************************************************/
  template<typename T>
  inline std::pair<bool, unsigned> Lib::SparseMatrix<T>::Slice::find(unsigned index) const
  {
    if (lastCheckedIndex == index) {
      return std::make_pair(data.size() > lastFoundIndex
                                && data[lastFoundIndex].first == index,
                            lastFoundIndex);
    }
    unsigned left = 0;
    unsigned right = data.size();
    while (left < right) {
      unsigned mid = left + (right - left) / 2;
      if (data[mid].first == index) {
        lastCheckedIndex = index;
        lastFoundIndex = mid;
        return std::make_pair(true, mid);
      }
      if (data[mid].first < index) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }
    lastCheckedIndex = index;
    lastFoundIndex = left;
    return std::make_pair(false, left);
  }

  template<typename T>
  inline void Lib::SparseMatrix<T>::Slice::set(unsigned index, T value)
  {
    auto [found, left] = find(index);
    if (found) {
      data[left].second = value;
      return;
    }
    data.insert(data.begin() + left, std::make_pair(index, value));
  }

  template<typename T>
  inline T Lib::SparseMatrix<T>::Slice::get(unsigned index, T def) const
  {
    auto [found, left] = find(index);
    if (found) {
      return data[left].second;
    }
    return def;
  }

  template<typename T>
  inline void SparseMatrix<T>::Slice::del(unsigned index)
  {
    auto [found, left] = find(index);
    if (found) {
      data.erase(data.begin() + left);
    }
  }

  /***********************************************************************************************/
  /*                                        SPARSE MATRIX                                        */
  /***********************************************************************************************/
  template<typename T>
  inline SparseMatrix<T>::SparseMatrix(unsigned rows, unsigned cols)
  {
    this->rows = rows;
    this->cols = cols;
    data.resize(rows);
  }

  template<typename T>
  inline void SparseMatrix<T>::set(unsigned row, unsigned col, T value)
  {
    data[row].set(col, value);
  }

  template<typename T>
  inline void SparseMatrix<T>::del(unsigned row, unsigned col)
  {
    data[row].del(col);
  }

  template<typename T>
  inline T SparseMatrix<T>::get(unsigned row, unsigned col, T def) const
  {
    return data[row].get(col, def);
  }

  template<typename T>
  inline std::vector<std::pair<unsigned, T>>& SparseMatrix<T>::getSetOnRow(unsigned row)
  {
    return data[row].data;
  }

  template<typename T>
  inline void SparseMatrix<T>::swapRows(unsigned row1, unsigned row2)
  {
    if (row1 == row2)
      return;
    std::swap(data[row1], data[row2]);
  }

  template<typename T>
  inline void SparseMatrix<T>::reset()
  {
    // keep the same capacity
    for (auto& row : data) {
      row.data.clear();
    }
  }

  template<typename T>
  inline void SparseMatrix<T>::reshape(unsigned rows, unsigned cols)
  {
    this->rows = rows;
    this->cols = cols;
    data.resize(rows);
  }
}
