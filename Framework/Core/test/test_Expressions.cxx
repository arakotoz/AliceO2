// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#define BOOST_TEST_MODULE Test Framework Expressions
#define BOOST_TEST_MAIN
#define BOOST_TEST_DYN_LINK

#include "Framework/Configurable.h"
#include "Framework/ExpressionHelpers.h"
#include "Framework/AnalysisDataModel.h"
#include "Framework/AODReaderHelpers.h"
#include <boost/test/unit_test.hpp>
#include <arrow/util/config.h>

using namespace o2::framework;
using namespace o2::framework::expressions;

namespace nodes
{
static BindingNode pt{"pt", 1, atype::FLOAT};
static BindingNode phi{"phi", 2, atype::FLOAT};
static BindingNode eta{"eta", 3, atype::FLOAT};

static BindingNode tgl{"tgl", 4, atype::FLOAT};
static BindingNode signed1Pt{"signed1Pt", 5, atype::FLOAT};
static BindingNode testInt{"testInt", 6, atype::INT32};
} // namespace nodes

namespace o2::aod::track
{
DECLARE_SOA_EXPRESSION_COLUMN(Pze, pz, float, o2::aod::track::tgl*(1.f / o2::aod::track::signed1Pt));
} // namespace o2::aod::track

BOOST_AUTO_TEST_CASE(TestTreeParsing)
{
  expressions::Filter f = ((nodes::phi > 1) && (nodes::phi < 2)) && (nodes::eta < 1);
  auto specs = createOperations(f);
  BOOST_REQUIRE_EQUAL(specs[0].left, (DatumSpec{1u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(specs[0].right, (DatumSpec{2u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(specs[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(specs[1].left, (DatumSpec{std::string{"eta"}, 3, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(specs[1].right, (DatumSpec{LiteralNode::var_t{1}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(specs[1].result, (DatumSpec{2u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(specs[2].left, (DatumSpec{3u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(specs[2].right, (DatumSpec{4u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(specs[2].result, (DatumSpec{1u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(specs[3].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(specs[3].right, (DatumSpec{LiteralNode::var_t{2}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(specs[3].result, (DatumSpec{4u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(specs[4].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(specs[4].right, (DatumSpec{LiteralNode::var_t{1}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(specs[4].result, (DatumSpec{3u, atype::BOOL}));

  expressions::Filter g = ((nodes::eta + 2.f) > 0.5) || ((nodes::phi - M_PI) < 3);
  auto gspecs = createOperations(g);
  BOOST_REQUIRE_EQUAL(gspecs[0].left, (DatumSpec{1u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(gspecs[0].right, (DatumSpec{2u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(gspecs[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(gspecs[1].left, (DatumSpec{3u, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(gspecs[1].right, (DatumSpec{LiteralNode::var_t{3}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(gspecs[1].result, (DatumSpec{2u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(gspecs[2].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(gspecs[2].right, (DatumSpec{LiteralNode::var_t{M_PI}, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(gspecs[2].result, (DatumSpec{3u, atype::DOUBLE}));

  BOOST_REQUIRE_EQUAL(gspecs[3].left, (DatumSpec{4u, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(gspecs[3].right, (DatumSpec{LiteralNode::var_t{0.5}, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(gspecs[3].result, (DatumSpec{1u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(gspecs[4].left, (DatumSpec{std::string{"eta"}, 3, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(gspecs[4].right, (DatumSpec{LiteralNode::var_t{2.f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(gspecs[4].result, (DatumSpec{4u, atype::FLOAT}));

  expressions::Filter h = (nodes::phi == 0) || (nodes::phi == 3);
  auto hspecs = createOperations(h);

  BOOST_REQUIRE_EQUAL(hspecs[0].left, (DatumSpec{1u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(hspecs[0].right, (DatumSpec{2u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(hspecs[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(hspecs[1].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(hspecs[1].right, (DatumSpec{LiteralNode::var_t{3}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(hspecs[1].result, (DatumSpec{2u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(hspecs[2].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(hspecs[2].right, (DatumSpec{LiteralNode::var_t{0}, atype::INT32}));
  BOOST_REQUIRE_EQUAL(hspecs[2].result, (DatumSpec{1u, atype::BOOL}));

  expressions::Filter u = nabs(nodes::eta) < 1.0 && nexp(nodes::phi + 2.0 * M_PI) > 3.0;
  auto uspecs = createOperations(std::move(u));
  BOOST_REQUIRE_EQUAL(uspecs[0].left, (DatumSpec{1u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(uspecs[0].right, (DatumSpec{2u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(uspecs[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(uspecs[1].left, (DatumSpec{3u, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(uspecs[1].right, (DatumSpec{LiteralNode::var_t{3.0}, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(uspecs[1].result, (DatumSpec{2u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(uspecs[2].left, (DatumSpec{4u, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(uspecs[2].right, (DatumSpec{}));
  BOOST_REQUIRE_EQUAL(uspecs[2].result, (DatumSpec{3u, atype::DOUBLE}));

  BOOST_REQUIRE_EQUAL(uspecs[3].left, (DatumSpec{std::string{"phi"}, 2, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(uspecs[3].right, (DatumSpec{LiteralNode::var_t{2.0 * M_PI}, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(uspecs[3].result, (DatumSpec{4u, atype::DOUBLE}));

  BOOST_REQUIRE_EQUAL(uspecs[4].left, (DatumSpec{5u, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(uspecs[4].right, (DatumSpec{LiteralNode::var_t{1.0}, atype::DOUBLE}));
  BOOST_REQUIRE_EQUAL(uspecs[4].result, (DatumSpec{1u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(uspecs[5].left, (DatumSpec{std::string{"eta"}, 3, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(uspecs[5].right, (DatumSpec{}));
  BOOST_REQUIRE_EQUAL(uspecs[5].result, (DatumSpec{5u, atype::FLOAT}));

  Configurable<float> pTCut{"pTCut", 0.5f, "Lower pT limit"};
  Filter ptfilter = o2::aod::track::pt > pTCut;
  BOOST_REQUIRE_EQUAL(ptfilter.node->self.index(), 2);
  BOOST_REQUIRE_EQUAL(ptfilter.node->left->self.index(), 1);
  BOOST_REQUIRE_EQUAL(ptfilter.node->right->self.index(), 3);
  auto ptfilterspecs = createOperations(ptfilter);
  BOOST_REQUIRE_EQUAL(ptfilterspecs[0].left, (DatumSpec{std::string{"fPt"}, typeid(o2::aod::track::Pt).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(ptfilterspecs[0].right, (DatumSpec{LiteralNode::var_t{0.5f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(ptfilterspecs[0].result, (DatumSpec{0u, atype::BOOL}));
}

BOOST_AUTO_TEST_CASE(TestGandivaTreeCreation)
{
  Projector pze = o2::aod::track::Pze::Projector();
  auto pzspecs = createOperations(pze);
  BOOST_REQUIRE_EQUAL(pzspecs[0].left, (DatumSpec{std::string{"fTgl"}, typeid(o2::aod::track::Tgl).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(pzspecs[0].right, (DatumSpec{1u, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(pzspecs[0].result, (DatumSpec{0u, atype::FLOAT}));

  BOOST_REQUIRE_EQUAL(pzspecs[1].left, (DatumSpec{LiteralNode::var_t{1.f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(pzspecs[1].right, (DatumSpec{std::string{"fSigned1Pt"}, typeid(o2::aod::track::Signed1Pt).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(pzspecs[1].result, (DatumSpec{1u, atype::FLOAT}));
  auto infield1 = o2::aod::track::Signed1Pt::asArrowField();
  auto infield2 = o2::aod::track::Tgl::asArrowField();
  auto resfield = o2::aod::track::Pze::asArrowField();
  auto schema = std::make_shared<arrow::Schema>(std::vector{infield1, infield2, resfield});
  auto gandiva_tree = createExpressionTree(pzspecs, schema);

  auto gandiva_expression = makeExpression(gandiva_tree, resfield);
  BOOST_CHECK_EQUAL(gandiva_expression->ToString(), "float multiply((float) fTgl, float divide((const float) 1 raw(3f800000), (float) fSigned1Pt))");
  auto projector = createProjector(schema, pzspecs, resfield);

  Projector pte = o2::aod::track::Pt::Projector();
  auto ptespecs = createOperations(pte);
  BOOST_REQUIRE_EQUAL(ptespecs[0].left, (DatumSpec{1u, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(ptespecs[0].right, (DatumSpec{}));
  BOOST_REQUIRE_EQUAL(ptespecs[0].result, (DatumSpec{0u, atype::FLOAT}));

  BOOST_REQUIRE_EQUAL(ptespecs[1].left, (DatumSpec{LiteralNode::var_t{1.f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(ptespecs[1].right, (DatumSpec{std::string{"fSigned1Pt"}, typeid(o2::aod::track::Signed1Pt).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(ptespecs[1].result, (DatumSpec{1u, atype::FLOAT}));

  auto infield3 = o2::aod::track::Signed1Pt::asArrowField();
  auto resfield2 = o2::aod::track::Pt::asArrowField();
  auto schema2 = std::make_shared<arrow::Schema>(std::vector{infield3, resfield2});
  auto gandiva_tree2 = createExpressionTree(ptespecs, schema2);

  auto gandiva_expression2 = makeExpression(gandiva_tree2, resfield2);
  BOOST_CHECK_EQUAL(gandiva_expression2->ToString(), "float absf(float divide((const float) 1 raw(3f800000), (float) fSigned1Pt))");

  auto projector_b = createProjector(schema2, ptespecs, resfield2);
  auto schema_p = o2::soa::createSchemaFromColumns(o2::aod::Tracks::persistent_columns_t{});
  auto projector_alt = o2::framework::expressions::createProjectors(o2::framework::pack<o2::aod::track::Pt>{}, schema_p);

  Filter bitwiseFilter = (o2::aod::track::flags & static_cast<uint32_t>(o2::aod::track::TPCrefit)) != 0u;
  auto bwf = createOperations(bitwiseFilter);
  BOOST_REQUIRE_EQUAL(bwf[0].left, (DatumSpec{1u, atype::UINT32}));
  BOOST_REQUIRE_EQUAL(bwf[0].right, (DatumSpec{LiteralNode::var_t{0u}, atype::UINT32}));
  BOOST_REQUIRE_EQUAL(bwf[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(bwf[1].left, (DatumSpec{std::string{"fFlags"}, typeid(o2::aod::track::Flags).hash_code(), atype::UINT32}));
  BOOST_REQUIRE_EQUAL(bwf[1].right, (DatumSpec{LiteralNode::var_t{static_cast<uint32_t>(o2::aod::track::TPCrefit)}, atype::UINT32}));
  BOOST_REQUIRE_EQUAL(bwf[1].result, (DatumSpec{1u, atype::UINT32}));

  auto infield4 = o2::aod::track::Flags::asArrowField();
  auto resfield3 = std::make_shared<arrow::Field>("out", arrow::boolean());
  auto schema_b = std::make_shared<arrow::Schema>(std::vector{infield4, resfield3});
  auto gandiva_tree3 = createExpressionTree(bwf, schema_b);
  BOOST_CHECK_EQUAL(gandiva_tree3->ToString(), "bool not_equal(uint32 bitwise_and((uint32) fFlags, (const uint32) 2), (const uint32) 0)");
  auto condition = expressions::makeCondition(gandiva_tree3);
  std::shared_ptr<gandiva::Filter> flt;
  auto s = gandiva::Filter::Make(schema_b, condition, &flt);
  BOOST_REQUIRE(s.ok());
}

BOOST_AUTO_TEST_CASE(TestConditionalExpressions)
{
  // simple conditional
  Filter cf = nabs(o2::aod::track::eta) < 1.0f && ifnode((o2::aod::track::pt < 1.0f), (o2::aod::track::phiraw > (float)(M_PI / 2.)), (o2::aod::track::phiraw < (float)(M_PI / 2.)));
  auto cfspecs = createOperations(cf);
  BOOST_REQUIRE_EQUAL(cfspecs[0].left, (DatumSpec{1u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(cfspecs[0].right, (DatumSpec{2u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(cfspecs[0].result, (DatumSpec{0u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[1].left, (DatumSpec{3u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(cfspecs[1].right, (DatumSpec{4u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(cfspecs[1].condition, (DatumSpec{5u, atype::BOOL}));
  BOOST_REQUIRE_EQUAL(cfspecs[1].result, (DatumSpec{2u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[2].left, (DatumSpec{std::string{"fPt"}, typeid(o2::aod::track::Pt).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[2].right, (DatumSpec{LiteralNode::var_t{1.0f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[2].result, (DatumSpec{5u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[3].left, (DatumSpec{std::string{"fRawPhi"}, typeid(o2::aod::track::RawPhi).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[3].right, (DatumSpec{LiteralNode::var_t{(float)(M_PI / 2.)}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[3].result, (DatumSpec{4u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[4].left, (DatumSpec{std::string{"fRawPhi"}, typeid(o2::aod::track::RawPhi).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[4].right, (DatumSpec{LiteralNode::var_t{(float)(M_PI / 2.)}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[4].result, (DatumSpec{3u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[5].left, (DatumSpec{6u, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[5].right, (DatumSpec{LiteralNode::var_t{1.0f}, atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[5].result, (DatumSpec{1u, atype::BOOL}));

  BOOST_REQUIRE_EQUAL(cfspecs[6].left, (DatumSpec{std::string{"fEta"}, typeid(o2::aod::track::Eta).hash_code(), atype::FLOAT}));
  BOOST_REQUIRE_EQUAL(cfspecs[6].right, (DatumSpec{}));
  BOOST_REQUIRE_EQUAL(cfspecs[6].result, (DatumSpec{6u, atype::FLOAT}));

  auto infield1 = o2::aod::track::Pt::asArrowField();
  auto infield2 = o2::aod::track::Eta::asArrowField();
  auto infield3 = o2::aod::track::RawPhi::asArrowField();
  auto schema = std::make_shared<arrow::Schema>(std::vector{infield1, infield2, infield3});
  auto gandiva_tree = createExpressionTree(cfspecs, schema);
  auto gandiva_condition = makeCondition(gandiva_tree);
  auto gandiva_filter = createFilter(schema, gandiva_condition);

  BOOST_CHECK_EQUAL(gandiva_tree->ToString(), "bool less_than(float absf((float) fEta), (const float) 1 raw(3f800000)) && if (bool less_than((float) fPt, (const float) 1 raw(3f800000))) { bool greater_than((float) fRawPhi, (const float) 1.5708 raw(3fc90fdb)) } else { bool less_than((float) fRawPhi, (const float) 1.5708 raw(3fc90fdb)) }");

  // nested conditional
  Filter cfn = o2::aod::track::signed1Pt > 0.f && ifnode(std::move(*cf.node), nabs(o2::aod::track::x) > 1.0f, nabs(o2::aod::track::y) > 1.0f);
  auto cfnspecs = createOperations(cfn);
  auto infield4 = o2::aod::track::Signed1Pt::asArrowField();
  auto infield5 = o2::aod::track::X::asArrowField();
  auto infield6 = o2::aod::track::Y::asArrowField();
  auto schema2 = std::make_shared<arrow::Schema>(std::vector{infield1, infield2, infield3, infield4, infield5, infield6});
  auto gandiva_tree2 = createExpressionTree(cfnspecs, schema2);
  auto gandiva_condition2 = makeCondition(gandiva_tree2);
  auto gandiva_filter2 = createFilter(schema2, gandiva_condition2);
  BOOST_REQUIRE_EQUAL(gandiva_tree2->ToString(),
                      "bool greater_than((float) fSigned1Pt, (const float) 0 raw(0)) && if (bool less_than(float absf((float) fEta), (const float) 1 raw(3f800000)) && if (bool less_than((float) fPt, (const float) 1 raw(3f800000))) { bool greater_than((float) fRawPhi, (const float) 1.5708 raw(3fc90fdb)) } else { bool less_than((float) fRawPhi, (const float) 1.5708 raw(3fc90fdb)) }) { bool greater_than(float absf((float) fX), (const float) 1 raw(3f800000)) } else { bool greater_than(float absf((float) fY), (const float) 1 raw(3f800000)) }");
}
