/*
 * Copyright (c) 2015-2021, The OSKAR Developers.
 * See the LICENSE file at the top-level directory of this distribution.
 */

#include <gtest/gtest.h>
#include "settings/oskar_settings_types.h"

using namespace oskar;

TEST(settings_types, Int)
{
    Int i;
    ASSERT_TRUE(i.set_default("1"));
    ASSERT_STREQ("1", i.get_default());
    ASSERT_STREQ("1", i.get_value());
    ASSERT_TRUE(i.is_default());
    ASSERT_TRUE(i.set_value("9"));
    ASSERT_FALSE(i.is_default());
    ASSERT_EQ(9, i.value());
    ASSERT_EQ(1, i.default_value());
    ASSERT_FALSE(i.set_value("2.111"));
    ASSERT_FALSE(i.set_value("foo"));

    // Comparison.
    {
        Int i1, i2;
        ASSERT_TRUE(i1.set_value("1"));
        ASSERT_TRUE(i2.set_value("1"));
        ASSERT_TRUE(i1 == i2);
        ASSERT_FALSE(i2 > i1);
        ASSERT_TRUE(i2.set_value("2"));
        ASSERT_FALSE(i1 == i2);
        ASSERT_TRUE(i2 > i1);
    }
}
