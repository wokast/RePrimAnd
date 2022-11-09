#include "datastore.h"
#include "hdf5store.h"
#include "interpol.h"
#include "star_sequence.h"
#include "tov_seqs_impl.h"
#include "star_seq_file.h"

namespace EOS_Toolkit {
namespace detail {

void star_seq_impl::save(datasink g) const
{
  g["mg_gm1"] = mg_gm1 * u.mass();
  g["mb_gm1"] = mb_gm1 * u.mass();
  g["rc_gm1"] = rc_gm1 * u.length();
  g["mi_gm1"] = mi_gm1 * u.mom_inertia();
  g["lt_gm1"] = lt_gm1;
  
  g["range_gm1"] = rg_gm1;
}

void star_branch_impl::save(datasink g) const
{
  g["xg_mg"] = xg_mg.rescale_x(u.mass());  
  g["range_gm1"] = rg_gm1;
  g["reference_gm1"] = gm1_ref;
  g["includes_max"] = incl_max;
}

auto star_seq_impl::from_datasource(datasource g, units u) 
-> std::shared_ptr<star_seq_impl>
{
  interpolator mg_gm1_si = g["mg_gm1"];
  interpolator mb_gm1_si = g["mb_gm1"];
  interpolator rc_gm1_si = g["rc_gm1"];
  interpolator mi_gm1_si = g["mi_gm1"];
  interpolator lt_gm1 = g["lt_gm1"];
  
  auto mg_gm1    = mg_gm1_si / u.mass(); 
  auto mb_gm1    = mb_gm1_si / u.mass(); 
  auto rc_gm1    = rc_gm1_si / u.length(); 
  auto mi_gm1    = mi_gm1_si / u.mom_inertia(); 
  
  interval<real_t> rg_gm1 = g["range_gm1"];

  return std::make_shared<star_seq_impl>(mg_gm1, mb_gm1, rc_gm1, 
                                         mi_gm1, lt_gm1, rg_gm1, u);
}

auto star_branch_impl::from_datasource(datasource g, units u) 
-> std::shared_ptr<star_branch_impl>
{
  interpolator xg_mg_si   = g["xg_mg"];  
  interval<real_t> rg_gm1 = g["range_gm1"];
  real_t gm1_ref          = g["reference_gm1"];
  bool incl_max           = g["includes_max"];  
  auto xg_mg              = xg_mg_si.rescale_x(1.0 / u.mass());
  
  return std::make_shared<star_branch_impl>(rg_gm1, xg_mg, 
                                                gm1_ref, incl_max, u);
  
}

auto load_star_seq(const datasource g, const units& u)
-> star_seq
{
  return star_seq(star_seq_impl::from_datasource(g, u));
}

auto load_star_branch(const datasource g, const units& u) 
-> star_branch
{
  auto sbr = star_branch_impl::from_datasource(g, u);
  auto seq = star_seq_impl::from_datasource(g / "star_sequence", u);

  return star_branch(seq, sbr);
}

}

void star_seq::save(datasink s) const
{
  implementation().save(s);  
}

void star_branch::save(datasink s) const
{
  implementation().save(s);  
  star_seq::implementation().save(s / "star_sequence");  
}

auto load_star_seq(std::string fname, const units& u)
-> star_seq
{
  auto g = make_hdf5_file_source(fname);
  return detail::load_star_seq(g / "star_sequence", u);
}

auto load_star_branch(std::string fname, const units& u)
-> star_branch
{
  auto g = make_hdf5_file_source(fname);
  return detail::load_star_branch(g / "star_sequence_branch", u);
}

void save_star_seq(std::string fname, const star_seq& seq)
{
  auto g = make_hdf5_file_sink(fname);
  seq.save(g / "star_sequence");
}

void save_star_branch(std::string fname, const star_branch& seq)
{
  auto g = make_hdf5_file_sink(fname);
  seq.save(g / "star_sequence_branch");  
}


} // namespace EOS_Toolkit
