#pragma once

#include <vector>

#include "Debug/Assertion.hpp"

#include <iostream>

namespace Lib {
  template <typename T>
  class SparseMatrix {
    private:
      class Slice {
        friend class SparseMatrix<T>;
        public:
          Slice() = default;

          void set(unsigned col, T value);

          /**
           * @brief Returns a pointer to the value at the given col. If no value is set at the given col, returns nullptr.
           * @note This function is not marked const because it modifies the lastCheckedIndex and lastFoundIndex members.
           */
          T* get(unsigned col);


          unsigned size() const;

        private:
          void del(unsigned col);

          void delCol(unsigned col);

          void clear();

          std::vector<std::pair<unsigned, T>> data;

          /**
           * @note This function is not marked const because it modifies the lastCheckedIndex and lastFoundIndex members.
           */
          std::pair<bool, unsigned> find(unsigned col) const;

          /**
           * @note Those fields are transparent to the user, so they are mutable.
           */
          mutable unsigned lastCheckedIndex = 0xFFFFFFFF;
          mutable unsigned lastFoundIndex = 0xFFFFFFFF;
      };
    public:
      SparseMatrix(unsigned rows, unsigned cols, T def=T());

      void set(unsigned row, unsigned col, T value);

      /**
       * @brief Returns a pointer to the value at the given col. If no value is set at the given col, returns nullptr.
       */
      T get(unsigned row, unsigned col) const;

      /**
       * @brief Returns a reference to the non-zero elements in the given row.
       */
      std::vector<std::pair<unsigned, T>>& getSetOnRow(unsigned row);

      void setDefault(T def);

      void swapRows(unsigned row1, unsigned row2);

      void delRow(unsigned row);

      void delCol(unsigned col);

      void clear();

      void reshape(unsigned rows, unsigned cols);

    private:
      /**
       * We cheat by making this mutable to keep the interface clean.
       * We need to do this because we use references to identify if an element was found, and the compiler does not see that we do not modify this field.
       */
      mutable std::vector<Slice> data;
      unsigned rows;
      unsigned cols;
      /**
       * default value
       */
      T defaultValue;

  };

  /***********************************************************************************************/
  /*                                            SLICE                                            */
  /***********************************************************************************************/
  template<typename T>
  inline std::pair<bool, unsigned> Lib::SparseMatrix<T>::Slice::find(unsigned col) const
  {
    if (lastCheckedIndex == col) {
      ASS_LE(lastFoundIndex, data.size());
      return std::make_pair(data.size() > lastFoundIndex
                                && data[lastFoundIndex].first == col,
                            lastFoundIndex);
    }
    unsigned left = 0;
    unsigned right = data.size();
    while (left < right) {
      unsigned mid = left + (right - left) / 2;
      if (data[mid].first == col) {
        lastCheckedIndex = col;
        lastFoundIndex = mid;
        return std::make_pair(true, mid);
      }
      if (data[mid].first < col) {
        left = mid + 1;
      } else {
        right = mid;
      }
    }
    lastCheckedIndex = col;
    lastFoundIndex = left;
    return std::make_pair(false, left);
  }

  template<typename T>
  inline void Lib::SparseMatrix<T>::Slice::set(unsigned col, T value)
  {
    auto [found, left] = find(col);
    if (found) {
      data[left].second = value;
      return;
    }
    data.insert(data.begin() + left, std::make_pair(col, value));
  }

  template<typename T>
  inline T* Lib::SparseMatrix<T>::Slice::get(unsigned col)
  {
    auto [found, left] = find(col);
    if (found)
      return &data[left].second;
    return nullptr;
  }

  template<typename T>
  inline void SparseMatrix<T>::Slice::del(unsigned col)
  {
    auto [found, left] = find(col);
    if (found)
      data.erase(data.begin() + left);
  }

  template<typename T>
  inline void SparseMatrix<T>::Slice::delCol(unsigned col)
  {
    auto [found, left] = find(col);
    if (found)
      data.erase(data.begin() + left);
    for (unsigned i = left; i < data.size(); i++)
      data[i].first--;
  }

  template<typename T>
  inline void SparseMatrix<T>::Slice::clear()
  {
    data.clear();
    lastCheckedIndex = 0xFFFFFFFF;
    lastFoundIndex = 0xFFFFFFFF;
  }

  /***********************************************************************************************/
  /*                                        SPARSE MATRIX                                        */
  /***********************************************************************************************/
  template<typename T>
  inline SparseMatrix<T>::SparseMatrix(unsigned rows, unsigned cols, T def)
  {
    this->rows = rows;
    this->cols = cols;
    this->defaultValue = def;

    data.resize(rows);
  }

  template<typename T>
  inline void SparseMatrix<T>::set(unsigned row, unsigned col, T value)
  {
    ASS_L(row, rows);
    ASS_L(col, cols);
    if (value == defaultValue)
      data[row].del(col);
    else
      data[row].set(col, value);
  }

  template<typename T>
  inline T SparseMatrix<T>::get(unsigned row, unsigned col) const
  {
    ASS_L(row, rows);
    ASS_L(col, cols);
    T* toReturn = data[row].get(col);
    return toReturn ? *toReturn : defaultValue;
  }

  template<typename T>
  inline std::vector<std::pair<unsigned, T>>& SparseMatrix<T>::getSetOnRow(unsigned row)
  {
    ASS_L(row, rows);
    return data[row].data;
  }

  template<typename T>
  inline void SparseMatrix<T>::setDefault(T def)
  {
    defaultValue = def;
  }

  template<typename T>
  inline void SparseMatrix<T>::swapRows(unsigned row1, unsigned row2)
  {
    ASS_L(row1, rows);
    ASS_L(row2, rows);
    if (row1 == row2)
      return;
    std::swap(data[row1], data[row2]);
  }

  template<typename T>
  inline void SparseMatrix<T>::delRow(unsigned row)
  {
    data.erase(data.begin() + row);
    rows--;
  }

  template<typename T>
  inline void SparseMatrix<T>::delCol(unsigned col)
  {
    for (auto& row : data) {
      row.delCol(col);
    }
    cols--;
  }

  template<typename T>
  inline void SparseMatrix<T>::clear()
  {
    // keep the same capacity
    for (auto& row : data) {
      row.clear();
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
