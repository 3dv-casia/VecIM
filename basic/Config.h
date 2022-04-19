#ifndef CONFIG_H
#define CONFIG_H

#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

namespace cm {

typedef boost::property_tree::ptree Config;

/*!
 * \brief Access the global singleton configuration.
 * Guaranteed to be destroyed and instantiated on first use.
 * https://stackoverflow.com/a/1008289/5487342
 * \return a boost property tree
 */
inline Config &get_config() {
  static Config config;
  return config;
}

/*!
 * \brief Read configurations from a xml file.
 * \param file a .xml file
 */
inline void read_config(const std::string &file) {
  Config &config = get_config();
  boost::property_tree::read_xml(file,
    config,
    boost::property_tree::xml_parser::trim_whitespace |
    boost::property_tree::xml_parser::no_comments);
}

/*!
 * \brief Write configurations to a xml file.
 * \param file a .xml file
 */
inline void write_config(const std::string &file) {
  const Config &config = get_config();
  boost::property_tree::write_xml(
    file,
    config,
    std::locale(),
    boost::property_tree::xml_writer_settings<std::string>(' ', 2));
}

} // namespace cm

#endif // CONFIG_H
