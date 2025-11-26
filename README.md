# plugin-SISYPHOS

Active GitHub repository for SISYPHOS, a tool for serial refinements in the Olex2 crystallographic software.

If you want to use SISYPHOS, you can download it as an extension module within the Olex2 crystallographic software suite.
Simply go to the HOME tab and click on "Extension modules", enter your details, and select the SISYPHPS module.

In order to run SISYPHOS in Olex2 you can copy this folder into the $OLEX2BASE\util\pyUtil\PluginLib\plugin-SISYPHOS (create the plugin-SISYPHOS folder there) and create a plugins.xld in the olex2 main folder which should have the following file content without the # symbols

'''
<Plugin
<SISYPHOS>
>
'''

You can also create your own plugins (use "spy.pt.make_new_plugin(NAME)" in the olex2 command line). For more info, visit [https://www.olexsys.org/olex2/docs/tasks/extending-olex2/plugin/]


