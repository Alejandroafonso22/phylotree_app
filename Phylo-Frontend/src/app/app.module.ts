import { NgModule } from '@angular/core';
import { BrowserModule } from '@angular/platform-browser';
import { AppRoutingModule } from './app-routing.module';
import { AppComponent } from './app.component';

// General imports
import { NgxPaginationModule } from 'ngx-pagination';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { ToastrModule } from 'ngx-toastr';
import { CookieService } from 'ngx-cookie-service';

// Layout imports
import { HomeComponent } from './layout/home/home.component';
import { FooterComponent } from './layout/footer/footer.component';
import { NavbarComponent } from './layout/navbar/navbar.component';
import { HeaderComponent } from './layout/header/header.component';
import { NotFoundComponent } from './meta/not-found/not-found.component';
import { FastaUploadComponent } from './layout/fasta-upload/fasta-upload.component';

// GOG component imports
import { GogHolderComponent } from './gog-holder/gog-holder.component';
import { MapViewComponent } from './gog-holder/map-view/map-view.component';
import { ShowUsrMarkersComponent } from './gog-holder/show-usr-markers/show-usr-markers.component';
import { AddEditUsrMarkersComponent } from './gog-holder/add-edit-usr-markers/add-edit-usr-markers.component';

// GOG services imports
import { MarkersService } from './services/gog/markers/markers.service';

// API consumption imports
import { HttpClientModule } from '@angular/common/http';
import { FormsModule, ReactiveFormsModule } from '@angular/forms';
import { DataVisualizationComponent } from './gog-holder/data-visualization/data-visualization.component';
import { MostRepeatedCountryComponent } from './gog-holder/data-visualization/most-repeated-country/most-repeated-country.component';
import { SpeciesWithMoreMarkersComponent } from './gog-holder/data-visualization/species-with-more-markers/species-with-more-markers.component';

// Validation Directives
// GOG directives
import { ValidateIntegerDirective } from './directives/gog/validate-integer/validate-integer.directive';
import { ValidateLongitudeDirective } from './directives/gog/validate-longitude/validate-longitude.directive';
import { ValidateLatitudeDirective } from './directives/gog/validate-latitude/validate-latitude.directive';
import { ValidateDateDirective } from './directives/gog/validate-date/validate-date.directive';
import { ValidateTimeDirective } from './directives/gog/validate-time/validate-time.directive';
import { ValidateCountryCodeDirective } from './directives/gog/validate-country-code/validate-country-code.directive';
import { ValidateStringDirective } from './directives/gog/validate-string/validate-string.directive';
import { DownloadNotificationsComponent } from './layout/navbar/download-notifications/download-notifications.component';

@NgModule({
  declarations: [
    AppComponent,
    GogHolderComponent,
    HomeComponent,
    FooterComponent,
    NavbarComponent,
    HeaderComponent,
    NotFoundComponent,
    MapViewComponent,
    ShowUsrMarkersComponent,
    AddEditUsrMarkersComponent,
    DataVisualizationComponent,
    MostRepeatedCountryComponent,
    SpeciesWithMoreMarkersComponent,
    FastaUploadComponent,
    ValidateIntegerDirective,
    ValidateLongitudeDirective,
    ValidateLatitudeDirective,
    ValidateDateDirective,
    ValidateTimeDirective,
    ValidateCountryCodeDirective,
    ValidateStringDirective,
    DownloadNotificationsComponent,
  ],
  imports: [
    BrowserModule,
    AppRoutingModule,
    HttpClientModule,
    FormsModule,
    ReactiveFormsModule,
    NgxPaginationModule,
    BrowserAnimationsModule,
    ToastrModule.forRoot()
  ],
  providers: [MarkersService, CookieService],
  bootstrap: [AppComponent]
})
export class AppModule { }
