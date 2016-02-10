/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

/** \file ptr_base.h
 * \brief Reference counted pointer to abstract base class for expressions.
 */

#ifndef SYMBOLIC_PTR_BASE_H

#define SYMBOLIC_PTR_BASE_H

#include <decl.h>
#include <algorithm>

namespace symbolic {
namespace internal {


/// Reference counted pointer to abstract base class for expressions.
class ptr_base {

  public:
    /// Constructor
    explicit ptr_base(ex_base *p);
    /// Copy constructor
    ptr_base(const ptr_base& p);
    /// Destructor
    ~ptr_base();

    /// Assignment
    ptr_base& operator=(const ptr_base& p);

    /// Reference
    ex_base& operator*() { return *m_p; }
    /// Const reference
    const ex_base& operator*() const { return *m_p; }
    /// Pointer.
    ex_base* operator->() { return m_p; }
    /// Const pointer.
    const ex_base* operator->() const { return m_p; }
    /// Pointer
    ex_base* get() { return m_p; }
    /// Const pointer.
    const ex_base* get() const { return m_p; }
    /// Pointer with static cast
    template <typename T> const T* get() const { return static_cast<const T*>(m_p); }

  private:
    // No default constructor
    ptr_base();
    ex_base *m_p;
    int *m_rc;
    void dec_rc();

}; /* class ptr_base */


} /* namespace internal */
} /* namespace symbolic */


#endif /* SYMBOLIC_PTR_BASE_H */
