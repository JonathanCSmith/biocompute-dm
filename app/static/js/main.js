/**
 * Created by jon on 04/11/15.
 */
function addNavigationBarBehaviours() {
    var navigationBar = $(".capture");

    /* affix the navbar after scroll below header */
    navigationBar.affix({
        offset: {
            top: $("header").height()
        }
    });

    /* Add the neccessary padding between nav and content to not skip content */
    navigationBar.on("affix.bs.affix", function() {
        $("main").css("paddingTop", "75px");
    });

    /* Remove the padding when it becomes unnecessary */
    navigationBar.on("affix-top.bs.affix", function() {
        $("main").css("padding-Top", "0");
    });
}

$(document).ready(function() {
    addNavigationBarBehaviours();
});
