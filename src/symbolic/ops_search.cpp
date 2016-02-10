/***********************************************************
 * (c) Kancelaria Prezesa Rady Ministrów 2012-2015         *
 * Treść licencji w pliku 'LICENCE'                        *
 *                                                         *
 * (c) Chancellery of the Prime Minister 2012-2015         *
 * License terms can be found in the file 'LICENCE'        *
 *                                                         *
 * Author: Grzegorz Klima                                  *
 ***********************************************************/

 /** \file ops.cpp
 * \brief Operations.
 */


#include <ex_num.h>
#include <ex_delta.h>
#include <ex_symb.h>
#include <ex_symbidx.h>
#include <ex_vart.h>
#include <ex_vartidx.h>
#include <ex_add.h>
#include <ex_mul.h>
#include <ex_pow.h>
#include <ex_func.h>
#include <ex_e.h>
#include <ex_sum.h>
#include <cmp.h>
#include <error.h>
#include <ops.h>
#include <iostream>
#include <cmath>
#include <climits>


using namespace symbolic;
using namespace symbolic::internal;




bool
symbolic::internal::has_Es(const ptr_base &p)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == DELTA) || (t == SYMB) || (t == SYMBIDX)
        || (t == VART) || (t == VARTIDX)) {
        return false;
    } else if (t == EX) {
        return true;
    } else if (t == ADD) {
        unsigned i, n = p.get<ex_add>()->no_ops();
        for (i = 0; i < n; ++i) {
            if (has_Es(p.get<ex_add>()->get_ops()[i].second))
                return true;
        }
        return false;
    } else if (t == MUL) {
        unsigned i, n = p.get<ex_mul>()->no_ops();
        for (i = 0; i < n; ++i) {
            if (has_Es(p.get<ex_mul>()->get_ops()[i].second))
                return true;
        }
        return false;
    } else if (t == POW) {
        const ex_pow *pt = p.get<ex_pow>();
        if (has_Es(pt->get_base())) return true;
        if (has_Es(pt->get_exp())) return true;
        return false;
    } else if (t == FUN) {
        return has_Es(p.get<ex_func>()->get_arg());
    } else if (t == SUM) {
        return has_Es(p.get<ex_sum>()->get_e());
    } else if (t == PROD) {
        return has_Es(p.get<ex_prod>()->get_e());
    } else if (t == IDX) {
        return has_Es(p.get<ex_idx>()->get_e());
    } else INTERNAL_ERROR
}



void
symbolic::internal::find_Es(const ptr_base &p, set_ex &sex)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == DELTA) || (t == SYMB) || (t == SYMBIDX)
        || (t == VART) || (t == VARTIDX)) {
        return;
    } else if (t == EX) {
        ptr_base pe = p.get<ex_e>()->get_arg();
        int l = p.get<ex_e>()->get_lag();
        sex.insert(ex(lag(pe, -l)));
        find_Es(pe, sex);
        return;
    } else if (t == ADD) {
        const num_ex_pair_vec &args = p.get<ex_add>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            find_Es(args[i].second, sex);
        }
        return;
    } else if (t == MUL) {
        const num_ex_pair_vec &args = p.get<ex_mul>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            find_Es(args[i].second, sex);
        }
        return;
    } else if (t == POW) {
        const ex_pow *pp = p.get<ex_pow>();
        find_Es(pp->get_base(), sex);
        find_Es(pp->get_exp(), sex);
        return;
    } else if (t == FUN) {
        find_Es(p.get<ex_func>()->get_arg(), sex);
        return;
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        set_ex s;
        find_Es(pt->get_e(), s);
        set_ex::const_iterator it;
        for (it = s.begin(); it != s.end(); ++it) {
            sex.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        set_ex s;
        find_Es(pt->get_e(), s);
        set_ex::const_iterator it;
        for (it = s.begin(); it != s.end(); ++it) {
            sex.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        set_ex s;
        find_Es(pt->get_e(), s);
        set_ex::const_iterator it;
        for (it = s.begin(); it != s.end(); ++it) {
            sex.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else INTERNAL_ERROR
}




void
symbolic::internal::collect(const ptr_base &p, set_ex &vars, set_ex &parms)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == DELTA)) {
        return;
    } else if (t == SYMB) {
        parms.insert(ex(p));
        return;
    } else if (t == SYMBIDX) {
        parms.insert(ex(p));
        return;
    } else if (t == VART) {
        vars.insert(ex(p.get<ex_vart>()->copy0()));
        return;
    } else if (t == VARTIDX) {
        vars.insert(ex(p.get<ex_vartidx>()->copy0()));
        return;
    } else if (t == EX) {
        ptr_base pe = p.get<ex_e>()->get_arg();
        collect(pe, vars, parms);
        return;
    } else if (t == ADD) {
        const num_ex_pair_vec &args = p.get<ex_add>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect(args[i].second, vars, parms);
        }
        return;
    } else if (t == MUL) {
        const num_ex_pair_vec &args = p.get<ex_mul>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect(args[i].second, vars, parms);
        }
        return;
    } else if (t == POW) {
        const ex_pow *pp = p.get<ex_pow>();
        collect(pp->get_base(), vars, parms);
        collect(pp->get_exp(), vars, parms);
        return;
    } else if (t == FUN) {
        collect(p.get<ex_func>()->get_arg(), vars, parms);
        return;
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        set_ex vs, ps;
        collect(pt->get_e(), vs, ps);
        set_ex::const_iterator it;
        for (it = vs.begin(); it != vs.end(); ++it) {
            vars.insert(ex(pt->get_ie(), *it));
        }
        for (it = ps.begin(); it != ps.end(); ++it) {
            parms.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        set_ex vs, ps;
        collect(pt->get_e(), vs, ps);
        set_ex::const_iterator it;
        for (it = vs.begin(); it != vs.end(); ++it) {
            vars.insert(ex(pt->get_ie(), *it));
        }
        for (it = ps.begin(); it != ps.end(); ++it) {
            parms.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        set_ex vs, ps;
        collect(pt->get_e(), vs, ps);
        set_ex::const_iterator it;
        for (it = vs.begin(); it != vs.end(); ++it) {
            vars.insert(ex(pt->get_ie(), *it));
        }
        for (it = ps.begin(); it != ps.end(); ++it) {
            parms.insert(ex(pt->get_ie(), *it));
        }
        return;
    } else INTERNAL_ERROR
}





void
symbolic::internal::collect_lags(const ptr_base &p, map_ex_int &map)
{
    if (!(p->flag () & HAST)) return;
    unsigned t = p->type();
    if (t == VART) {
        int l;
        l = p.get<ex_vart>()->get_lag();
        if (l >= -1) return;
        ex v = ex(p.get<ex_vart>()->copy0());
        map_ex_int::iterator it = map.find(v);
        if (it != map.end()) {
            if (it->second > l) it->second = l;
        } else {
            map.insert(std::pair<ex, int>(v, l));
        }
        return;
    } else if (t == VARTIDX) {
        int l;
        l = p.get<ex_vartidx>()->get_lag();
        if (l >= -1) return;
        ex v = ex(p.get<ex_vartidx>()->copy0());
        map_ex_int::iterator it = map.find(v);
        if (it != map.end()) {
            if (it->second > l) it->second = l;
        } else {
            map.insert(std::pair<ex, int>(v, l));
        }
        return;
    } else if (t == EX) {
        ptr_base pe = p.get<ex_e>()->get_arg();
        collect_lags(pe, map);
        return;
    } else if (t == ADD) {
        const num_ex_pair_vec &args = p.get<ex_add>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_lags(args[i].second, map);
        }
        return;
    } else if (t == MUL) {
        const num_ex_pair_vec &args = p.get<ex_mul>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_lags(args[i].second, map);
        }
        return;
    } else if (t == POW) {
        const ex_pow *pp = p.get<ex_pow>();
        collect_lags(pp->get_base(), map);
        collect_lags(pp->get_exp(), map);
        return;
    } else if (t == FUN) {
        collect_lags(p.get<ex_func>()->get_arg(), map);
        return;
    } else if (t == SUM) {
        map_ex_int temp;
        const ex_sum *pt = p.get<ex_sum>();
        collect_lags(pt->get_e(), temp);
        map_ex_int::const_iterator it = temp.begin(), itt;
        for (; it != temp.end(); ++it) {
            ex ins(ex_idx::create(pt->get_ie(), it->first.m_ptr));
            map_ex_int::iterator itt = map.find(ins);
            if (itt != map.end()) {
                if (itt->second > it->second) itt->second = it->second;
            } else {
                map.insert(std::pair<ex, int>(ins, it->second));
            }
        }
        return;
    } else if (t == PROD) {
        map_ex_int temp;
        const ex_prod *pt = p.get<ex_prod>();
        collect_lags(pt->get_e(), temp);
        map_ex_int::const_iterator it = temp.begin(), itt;
        for (; it != temp.end(); ++it) {
            ex ins(ex_idx::create(pt->get_ie(), it->first.m_ptr));
            map_ex_int::iterator itt = map.find(ins);
            if (itt != map.end()) {
                if (itt->second > it->second) itt->second = it->second;
            } else {
                map.insert(std::pair<ex, int>(ins, it->second));
            }
        }
        return;
    } else if (t == IDX) {
        map_ex_int temp;
        const ex_idx *pt = p.get<ex_idx>();
        collect_lags(pt->get_e(), temp);
        map_ex_int::const_iterator it = temp.begin(), itt;
        for (; it != temp.end(); ++it) {
            ex ins(ex_idx::create(pt->get_ie(), it->first.m_ptr));
            map_ex_int::iterator itt = map.find(ins);
            if (itt != map.end()) {
                if (itt->second > it->second) itt->second = it->second;
            } else {
                map.insert(std::pair<ex, int>(ins, it->second));
            }
        }
        return;
    } else INTERNAL_ERROR
}





void
symbolic::internal::collect_idx(const ptr_base &p, std::set<unsigned> &iset)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == SYMB) || (t == VART)) {
        return;
    } else if (t == DELTA) {
        const ex_delta *pt = p.get<ex_delta>();
        int i;
        i = pt->m_idx1;
        if (i > 0) iset.insert(i);
        i = pt->m_idx2;
        if (i > 0) iset.insert(i);
        return;
    } else if (t == SYMBIDX) {
        const ex_symbidx *pt = p.get<ex_symbidx>();
        int i;
        i = pt->m_idx1;
        if (i > 0) iset.insert(i);
        i = pt->m_idx2;
        if (i > 0) iset.insert(i);
        i = pt->m_idx3;
        if (i > 0) iset.insert(i);
        i = pt->m_idx4;
        if (i > 0) iset.insert(i);
        return;
    } else if (t == VARTIDX) {
        const ex_vartidx *pt = p.get<ex_vartidx>();
        int i;
        i = pt->m_idx1;
        if (i > 0) iset.insert(i);
        i = pt->m_idx2;
        if (i > 0) iset.insert(i);
        i = pt->m_idx3;
        if (i > 0) iset.insert(i);
        i = pt->m_idx4;
        if (i > 0) iset.insert(i);
        return;
    } else if (t == ADD) {
        const num_ex_pair_vec &args = p.get<ex_add>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_idx(args[i].second, iset);
        }
        return;
    } else if (t == MUL) {
        const num_ex_pair_vec &args = p.get<ex_mul>()->get_ops();
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_idx(args[i].second, iset);
        }
        return;
    } else if (t == POW) {
        collect_idx(p.get<ex_pow>()->get_base(), iset);
        collect_idx(p.get<ex_pow>()->get_exp(), iset);
        return;
    } else if (t == FUN) {
        collect_idx(p.get<ex_func>()->get_arg(), iset);
        return;
    } else if (t == EX) {
        collect_idx(p.get<ex_e>()->get_arg(), iset);
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        int i;
        i = pt->get_ie().get_excl_id();
        if (i > 0) iset.insert(i);
        std::set<unsigned> is;
        collect_idx(pt->get_e(), is);
        is.erase(pt->get_ie().get_id());
        for (std::set<unsigned>::const_iterator it = is.begin(); it != is.end(); ++it)
            iset.insert(*it);
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        int i;
        i = pt->get_ie().get_excl_id();
        if (i > 0) iset.insert(i);
        std::set<unsigned> is;
        collect_idx(pt->get_e(), is);
        is.erase(pt->get_ie().get_id());
        for (std::set<unsigned>::const_iterator it = is.begin(); it != is.end(); ++it)
            iset.insert(*it);
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        int i;
        i = pt->get_ie().get_excl_id();
        if (i > 0) iset.insert(i);
        std::set<unsigned> is;
        collect_idx(pt->get_e(), is);
        is.erase(pt->get_ie().get_id());
        for (std::set<unsigned>::const_iterator it = is.begin(); it != is.end(); ++it)
            iset.insert(*it);
    } else INTERNAL_ERROR
}



void
symbolic::internal::collect_fidx(const ptr_base &p, std::map<unsigned, unsigned> &im)
{
    unsigned t = p->type();
    if ((t == NUM) || (t == SYMB) || (t == VART)|| (t == DELTA)
        || (t == SYMBIDX) || (t == VARTIDX)) {
        return;
    } else if (t == ADD) {
        const num_ex_pair_vec &args = p.get<ex_add>()->get_ops();
        std::map<unsigned, unsigned> imcpy = im;
        std::map<unsigned, unsigned>::iterator it;
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_fidx(args[i].second, im);
            for (it = im.begin(); it != im.end(); ++it) {
                if ((it->second == 0) && (imcpy.find(it->first) == imcpy.end()))
                    im.erase(it);
            }
        }
        return;
    } else if (t == MUL) {
        const num_ex_pair_vec &args = p.get<ex_mul>()->get_ops();
        std::map<unsigned, unsigned> imcpy = im;
        std::map<unsigned, unsigned>::iterator it;
        unsigned i, n = args.size();
        for (i = 0; i < n; ++i) {
            collect_fidx(args[i].second, im);
            for (it = im.begin(); it != im.end(); ++it) {
                if ((it->second == 0) && (imcpy.find(it->first) == imcpy.end()))
                    im.erase(it);
            }
        }
        return;
    } else if (t == POW) {
        std::map<unsigned, unsigned> imcpy = im;
        std::map<unsigned, unsigned>::iterator it;
        collect_fidx(p.get<ex_pow>()->get_base(), im);
        for (it = im.begin(); it != im.end(); ++it) {
            if ((it->second == 0) && (imcpy.find(it->first) == imcpy.end()))
                im.erase(it);
        }
        collect_fidx(p.get<ex_pow>()->get_exp(), im);
        for (it = im.begin(); it != im.end(); ++it) {
            if ((it->second == 0) && (imcpy.find(it->first) == imcpy.end()))
                im.erase(it);
        }
        return;
    } else if (t == FUN) {
        collect_fidx(p.get<ex_func>()->get_arg(), im);
        return;
    } else if (t == EX) {
        collect_fidx(p.get<ex_e>()->get_arg(), im);
    } else if (t == SUM) {
        const ex_sum *pt = p.get<ex_sum>();
        symbolic::internal::collect_fidx(pt->get_ie(), im);
        symbolic::internal::collect_fidx(pt->get_e(), im);
    } else if (t == PROD) {
        const ex_prod *pt = p.get<ex_prod>();
        symbolic::internal::collect_fidx(pt->get_ie(), im);
        symbolic::internal::collect_fidx(pt->get_e(), im);
    } else if (t == IDX) {
        const ex_idx *pt = p.get<ex_idx>();
        symbolic::internal::collect_fidx(pt->get_ie(), im);
        symbolic::internal::collect_fidx(pt->get_e(), im);
    } else INTERNAL_ERROR
}


void
symbolic::internal::collect_fidx(const idx_ex &ie, std::map<unsigned, unsigned> &im)
{
    if (!ie) return;
    std::map<unsigned, unsigned>::iterator it;
    unsigned i = ie.get_id();
    if ((it = im.find(i)) != im.end()) { ++it->second; } else { im[i] = 0; }
}














